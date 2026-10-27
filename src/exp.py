# evndyn
import torch
from module import scVAE
from dataset import VaeDataSet, VaeDataManager
from torch.utils.data import DataLoader
import numpy as np
import copy
from utils import StageChange, DecreasingTemp, temp_scheduler
from torch.utils.tensorboard import SummaryWriter
from tqdm import tqdm
import sys
import utils
import os

class VaeExperiment:
    def __init__(self, adata, adata_mg, model_params, lr, lr_scheduler, val_ratio, test_ratio, batch_size, num_workers, loss_hparams, save_model_path, epoch_num, sample_info, tb_writer=None, warmup_epochs=3000, max_trainx_epochs=None,
                 stage_patience_x=10, stage_patience_pi=20, lr_patience=10, lr_factor=0.1):
        self.warmup_epochs = warmup_epochs
        self.edm = VaeDataManager(adata, adata_mg, batch_size, num_workers, val_ratio, test_ratio, sample_info=sample_info)
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        # train_pi/train_sp で enc_z(z) へ勾配を流すか。True=現行挙動(凍結), False=解凍。
        self.freeze_z_after_x = model_params.pop('freeze_z_after_x', True)
        # enc_z を凍結した train_pi/train_sp では q(z|x) が変わらないので、段階の切り替え時に一度だけ計算して保持し、
        # 以後は x を読まずに学習する (_build_qz_cache)。False なら旧来どおり毎バッチ x を読んで再構成項も計算する。
        self.cache_qz_after_x = model_params.pop('cache_qz_after_x', True)
        self._qz_cache = None
        self.model = scVAE(**model_params)
        self.model.to(self.device)
        # xnorm を毎バッチ dense 転送せず GPU で rank-1 再構成するための係数(1回だけ転送)
        self._col_mean_dev = self.edm.col_mean.to(self.device)
        self._scale = float(self.edm.scale)
        self.optimizer = torch.optim.Adam(self.model.parameters(), lr=lr)
        self.train_loss_list = []
        self.test_loss_list = []
        self.best_loss = None
        self.loss_hparams = loss_hparams
        self.turn_on_sp = loss_hparams.get('lmd_pi', 0) != 0 or loss_hparams.get('lmd_phi', 0) != 0
        self.save_model_path = save_model_path
        self.tb_writer = tb_writer
        self.epoch_num = epoch_num
        # 早期停止(stage遷移)とLRスケジューラの patience は workflow から制御可能。
        # stage_patience_pi は train_pi / train_sp 共通。LR patience/factor は ReduceLROnPlateau 用。
        self.stage_patience_x = stage_patience_x
        self.stage_patience_pi = stage_patience_pi
        self.lr_patience = lr_patience
        self.lr_factor = lr_factor
        self.stagemanager = StageChange(self.save_model_path, patience=self.stage_patience_x)
        if lr_scheduler == True:
            self.lr_scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(self.optimizer, 'min', factor=self.lr_factor, patience=self.lr_patience)
        else:
            self.lr_scheduler = None
        self.lr = lr
        print('model_params:', model_params)
        if max_trainx_epochs != None:
            self.max_trainx_epochs = max_trainx_epochs
        else:
            self.max_trainx_epochs = epoch_num

    def train_epoch(self, epoch, stage):
        self.model.train()
        total_loss = 0
        entry_num = 0
        
        # Apply warmup learning rate
        # 段階 1 (train_x) だけに適用する。旧版は段階 1 が warmup_epochs より早く終わると、段階 2 の開始時に
        # update_optimizer_lr で lr*0.1 に下げた値と ReduceLROnPlateau の引き下げを毎エポック上書きしていた (2026-10-10 修正)
        if stage == 'train_x' and epoch < self.warmup_epochs:
            warmup_factor = min(1.0, epoch / self.warmup_epochs)
            for param_group in self.optimizer.param_groups:
                param_group['lr'] = self.lr * warmup_factor
            
        _prof = os.environ.get('CLONEVI_PROFILE') == '1'
        if _prof:
            import time as _time
            _t_data = 0.0; _t_comp = 0.0; _t0 = _time.time()
        if self._qz_cache is not None:
            loss_list = self._train_batches_from_qz(stage)
        else:
            for x, row_mean_b, d, b, sample_info in self.edm.train_loader:
                if _prof:
                    _tm = _time.time(); _t_data += _tm - _t0

                x = x.to(self.device)
                row_mean_b = row_mean_b.to(self.device)
                # xnorm を GPU 側で rank-1 再構成。(scale*row_mean[i])*col_mean[j] と演算順序が一致し
                # 結果はビット単位で同一（旧 __getitem__ の dense 版と完全等価）。
                x_norm_mat = (self._scale * row_mean_b).unsqueeze(1) * self._col_mean_dev
                d = d.to(self.device)
                b = b.to(self.device)
                sample_info = sample_info.to(self.device)
                self.optimizer.zero_grad()
                loss_list = self.model.elbo_loss(
                    x, x_norm_mat, d, b, sample_info, self.loss_hparams)
                if stage == 'train_x':
                    normalized_loss = loss_list[0]
                elif stage == 'train_pi':
                    normalized_loss = loss_list[1]
                elif stage == 'train_sp':
                    normalized_loss = loss_list[2]
                normalized_loss.backward()
                self.optimizer.step()
                if _prof:
                    torch.cuda.synchronize(); _t0 = _time.time(); _t_comp += _t0 - _tm
        if _prof:
            print(f'[profile] data_wait={_t_data:.1f}s gpu_compute={_t_comp:.1f}s', flush=True)
        loss = loss_list[2]
        if self.tb_writer != None:
            self.tb_writer.add_scalar('all_loss/train', loss, epoch)
            self.tb_writer.add_scalar('kld/train', loss_list[3], epoch)
            if self._qz_cache is None: # キャッシュ使用中は x の再構成項を計算していない
                self.tb_writer.add_scalar('x/train', loss_list[4], epoch)
            self.tb_writer.add_scalar('pi/train', loss_list[5], epoch)
            self.tb_writer.add_scalar('sp_pi/train', loss_list[6], epoch)
            self.tb_writer.add_scalar('sp_phi/train', loss_list[7], epoch)
            self.tb_writer.add_scalars('Temperatures', {'kld_temp':self.loss_hparams['kld_temp'].temp,
                                                        'binom_temp':self.loss_hparams['binom_temp'].temp,
                                                        'sp_temp':self.loss_hparams['sp_temp'].temp}, epoch)
            #print('T_kld:', self.loss_hparams['T_kld'].temp)
            #print('T_binom:', self.loss_hparams['T_binom'].temp)
        
            
            
            # Log current learning rate
            current_lr = self.optimizer.param_groups[0]['lr']
            self.tb_writer.add_scalar("learning_rate", current_lr, epoch)
        return(loss_list)
        #return(total_loss / entry_num)
        #return(total_loss)

    def _eval_batched(self, X, Xn, D, B, SI, bs=512):
        # 評価をミニバッチ化(no_grad)。val/test全体(数千×85859)を一括forwardするとNB log_probの
        # 中間テンソルでGPU OOMするため。各loss成分はper-cell meanなのでバッチサイズ重み付き平均で
        # 全体平均に一致(監視用指標)。
        self.model.eval()
        torch.cuda.empty_cache()
        N = X.shape[0]
        acc = None; total = 0
        with torch.no_grad():
            for s0 in range(0, N, bs):
                s1 = min(s0 + bs, N); n = s1 - s0
                x = X[s0:s1].to(self.device)
                xn = Xn[s0:s1].to(self.device)
                d = D[s0:s1].to(self.device)
                b = B[s0:s1].to(self.device)
                si = SI[s0:s1].to(self.device)
                ll = [v.detach() for v in self.model.elbo_loss(x, xn, d, b, si, self.loss_hparams)]
                acc = [v * n for v in ll] if acc is None else [a + v * n for a, v in zip(acc, ll)]
                total += n
                del x, xn, d, b, si, ll
        return [a / total for a in acc]

    def _build_qz_cache(self):
        # enc_z を凍結した後の q(z|x) の平均・標準偏差を train/val の全細胞について一度だけ計算し、
        # d, b, sample_info と共に GPU に置く。encoder には dropout も batch 依存の正規化もないので、
        # 毎バッチ x から計算する値と同じ (行列積の丸め誤差を除く)。
        # x のモデル (enc_z, dec_z2x, theta) はここで確定するので logtheta_x も凍結する
        # (旧来の経路では train_pi 以降も x の再構成項を通じて logtheta_x だけが更新されていた。下流解析では未使用)
        self.model.eval()
        edm = self.edm
        def _encode(rows, bs=edm.train_loader.batch_size): # x を密行列にする行数は学習のバッチと同じ (メモリ上限を揃える)
            rows = np.asarray(rows)
            mus, scales = [], []
            with torch.no_grad():
                for s0 in range(0, len(rows), bs):
                    r = rows[s0:s0 + bs]
                    x = torch.from_numpy(np.asarray(edm.x[r].todense()).astype(np.float32)).to(self.device)
                    mu, scale = self.model.encode_qz(x, edm.sample_info[r].to(self.device))
                    mus.append(mu); scales.append(scale)
            return (torch.cat(mus), torch.cat(scales))
        self._store_qz_cache(_encode(edm.train_idx), _encode(edm.val_idx))

    def set_qz_cache(self, qz_mu, qz_scale=None):
        # 段階 1 を学習せずに段階 2 から始めるとき (workflow.estimate_probability の qz 引数)。
        # 全細胞 (adata の行順) の q(z|x) の平均・標準偏差を受け取る。qz_scale=None なら z を平均に固定する (外部の埋め込みなど)。
        # z が与えられるので encoder は使わない (凍結扱い)
        self.freeze_z_after_x = True
        mu = torch.as_tensor(np.asarray(qz_mu), dtype=torch.float32)
        scale = None if qz_scale is None else torch.as_tensor(np.asarray(qz_scale), dtype=torch.float32)
        def _rows(rows):
            rows = torch.as_tensor(np.asarray(rows))
            return (mu[rows].to(self.device), None if scale is None else scale[rows].to(self.device))
        self._store_qz_cache(_rows(self.edm.train_idx), _rows(self.edm.val_idx))

    def _store_qz_cache(self, train_qz, val_qz):
        # train/val の (平均, 標準偏差) を d, b, sample_info と共に GPU に置く
        edm = self.edm
        self.model.logtheta_x.requires_grad = False
        train_idx = np.asarray(edm.train_idx)
        self._qz_cache = {
            'train': (*train_qz, edm.d[train_idx].to(self.device), edm.b[train_idx].to(self.device), edm.sample_info[train_idx].to(self.device)),
            'val': (*val_qz, edm.val_d.to(self.device), edm.val_b.to(self.device), edm.val_sample_info.to(self.device))}
        # 添字だけを回す DataLoader。train_loader と同じ RandomSampler なので毎エポックの乱数の引き方
        # (base_seed と sampler の seed) が同じで、ミニバッチの並びも x を読む場合と一致する
        self._qz_idx_loader = DataLoader(range(len(train_idx)), batch_size=edm.train_loader.batch_size, shuffle=True)
        print(f'cached q(z|x): train {len(train_idx)} cells, val {len(edm.val_idx)} cells')

    def _train_batches_from_qz(self, stage):
        # train_epoch の x を読まない版 (enc_z 凍結後の train_pi / train_sp)
        assert stage in ('train_pi', 'train_sp')
        mu, scale, d, b, si = self._qz_cache['train']
        for idx in self._qz_idx_loader:
            idx = idx.to(self.device)
            self.optimizer.zero_grad()
            loss_list = self.model.elbo_loss_from_qz(mu[idx], None if scale is None else scale[idx], d[idx], b[idx], si[idx], self.loss_hparams)
            normalized_loss = loss_list[1] if stage == 'train_pi' else loss_list[2]
            normalized_loss.backward()
            self.optimizer.step()
        return(loss_list)

    def _eval_from_qz(self, MU, SC, D, B, SI, bs=512):
        # _eval_batched の x を読まない版。z を引く乱数の消費を揃えるため同じ 512 細胞ずつに分ける
        self.model.eval()
        N = MU.shape[0]
        acc = None; total = 0
        with torch.no_grad():
            for s0 in range(0, N, bs):
                s1 = min(s0 + bs, N); n = s1 - s0
                ll = [v.detach() for v in self.model.elbo_loss_from_qz(MU[s0:s1], None if SC is None else SC[s0:s1], D[s0:s1], B[s0:s1], SI[s0:s1], self.loss_hparams)]
                acc = [v * n for v in ll] if acc is None else [a + v * n for a, v in zip(acc, ll)]
                total += n
        return [a / total for a in acc]

    def evaluate(self):
        if self._qz_cache is not None:
            return self._eval_from_qz(*self._qz_cache['val'])
        return self._eval_batched(self.edm.val_x, self.edm.val_xnorm_mat,
                                  self.edm.val_d, self.edm.val_b, self.edm.val_sample_info)

    def test(self):
        loss_list = self._eval_batched(self.edm.test_x, self.edm.test_xnorm_mat,
                                       self.edm.test_d, self.edm.test_b, self.edm.test_sample_info)
        return loss_list[2] / self.edm.test_x.shape[0]

    def train_total(self, train_mode):
        if train_mode == 'train_all' or train_mode == 'train_x':
            stage = 'train_x'
            for param in self.model.dec_z2pi.parameters(): # z→piのパラメータを固定
                param.requires_grad = False
        elif train_mode == 'train_pi':
            stage = 'train_pi'
            for param in self.model.dec_z2pi.parameters(): # z→piのパラメータを解凍
                param.requires_grad = True
            for param in self.model.dec_z2x.parameters(): # z→xのパラメータを凍結
                param.requires_grad = False
            for param in self.model.enc_z.parameters(): # zを固定 (freeze_z_after_x=Falseで解凍し勾配を流す)
                param.requires_grad = not self.freeze_z_after_x
            if self.freeze_z_after_x and self.cache_qz_after_x and self._qz_cache is None: # 以後 x を読まない (set_qz_cache 済みならそれを使う)
                self._build_qz_cache()
            # 段階 1 から移るときと同じ条件で始める (学習率 × 0.1、新しい Adam、段階 2 の patience)
            self.init_optimizer(self.lr * 0.1)
            self.stagemanager = StageChange(self.save_model_path, patience=self.stage_patience_pi)
            self.update_optimizer_lr(self.lr * 0.1)
        elif train_mode == 'train_sp':
            stage = 'train_sp'
        else:
            print('invalid train mode.')
            sys.exit(1)
        print(f'start {stage}')


        self.completed = True # 最後の段階が早期停止で終わったか (epoch_num 到達・NaN なら False)
        self.stage_log = []
        for epoch in tqdm(range(self.epoch_num), leave=False, ascii=True):
            state_dict = copy.deepcopy(self.model.state_dict())
            loss_list = self.train_epoch(epoch, stage)
            loss = loss_list[2].item()
            #print(loss)
            self.train_loss_list.append(loss)
            if np.isnan(loss):
                self.model.load_state_dict(state_dict)
                print(f'WARNING: loss is NaN at epoch {epoch} ({stage}); training stopped')
                self.completed = False
                break

            #if epoch % 2000 == 0:
            #    self._sample_and_save_z(epoch, stage)
            #    self.evaluate_z_clustering(epoch, stage)

            if epoch % 10 == 0:
                val_loss_list = self.evaluate()
                val_loss = val_loss_list[2].item()
                self.test_loss_list.append(val_loss)
                if self.tb_writer != None:
                    self.tb_writer.add_scalar('all_loss/val', val_loss, epoch)
                    self.tb_writer.add_scalar('kld/val', val_loss_list[3], epoch)
                    if self._qz_cache is None: # キャッシュ使用中は x の再構成項を計算していない
                        self.tb_writer.add_scalar('x/val', val_loss_list[4], epoch)
                    self.tb_writer.add_scalar('pi/val', val_loss_list[5], epoch)
                    self.tb_writer.add_scalar('sp_pi/val', val_loss_list[6], epoch)
                    self.tb_writer.add_scalar('sp_phi/val', val_loss_list[7], epoch)
                #print(f'loss at epoch {epoch} is {loss}')
                #print(f'val loss at epoch {epoch} is {val_loss}')
                #stagemanager(val_loss_list, self.model, stage)
                # 最良の重みは段階のすべての評価から選ぶ (温度で重みを付けない検証損失、StageChange.criterion)。
                # 段階を終える判定と学習率の引き下げは、その段階の warm-up と温度の下降が終わってから (_schedule_done)
                ready = self._schedule_done(epoch, stage)
                self.stagemanager(val_loss_list, self.model, stage, epoch=epoch, can_stop=ready)
                if stage == 'train_x':
                    if self.lr_scheduler != None and ready:
                        self.lr_scheduler.step(StageChange.criterion(val_loss_list, stage))
                    if self.stagemanager.flag==True or epoch>=self.max_trainx_epochs:
                        if train_mode == 'train_all':
                            self._restore_best(stage, epoch) # 段階 1 で検証損失が最良だった重みから段階 2 を始める
                            self.init_optimizer(self.lr * 0.1) # 戻した重みに後半のモーメントを持ち込まないよう Adam を作り直す
                            self.stagemanager = StageChange(self.save_model_path, patience=self.stage_patience_pi)
                            self.update_optimizer_lr(self.lr * 0.1) # 学習率スケジューラも新しい optimizer に付け直す
                            stage = 'train_pi'
                            for param in self.model.dec_z2pi.parameters(): # z→piのパラメータを解凍
                                param.requires_grad = True
                            for param in self.model.dec_z2x.parameters(): # z→xのパラメータを凍結
                                param.requires_grad = False
                            for param in self.model.enc_z.parameters(): # zを固定 (freeze_z_after_x=Falseで解凍し勾配を流す)
                                param.requires_grad = not self.freeze_z_after_x
                            if self.freeze_z_after_x and self.cache_qz_after_x: # 以後 x を読まない
                                self._build_qz_cache()
                            print(f'start {stage}')
                        else:
                            print('end training')
                            break
                elif stage == 'train_pi':
                    if self.lr_scheduler != None and ready:
                        self.lr_scheduler.step(val_loss_list[5])
                    if self.stagemanager.flag==True:
                        if self.turn_on_sp==True:
                            self._restore_best(stage, epoch) # 段階 2 で検証損失が最良だった重みから段階 3 を始める
                            self.init_optimizer(self.lr * 0.1) # 戻した重みに後半のモーメントを持ち込まないよう Adam を作り直す
                            self.stagemanager = StageChange(self.save_model_path, patience=self.stage_patience_pi)
                            stage = 'train_sp'
                            self.update_optimizer_lr(self.lr * 0.1)
                            print(f'start {stage}')
                        else:
                            print('end training')
                            break
                elif stage == 'train_sp':
                    if self.lr_scheduler != None and ready:
                        self.lr_scheduler.step(val_loss_list[5]+val_loss_list[6]+val_loss_list[7])
                    if self.stagemanager.flag==True:
                        print('end training')
                        break
                else:
                    print('invalid stage name.')
                    sys.exit(1)

                # Learning rate logging is now done in train_epoch method
            ''''  
                if stage == 'train_sp':
                    temp_scheduler(self.loss_hparams, stage, self.turn_on_sp, val_loss_list[5])
            
            if stage == 'train_x' or stage == 'train_pi':

                temp_scheduler(self.loss_hparams, stage, self.turn_on_sp)
            '''
            temp_scheduler(self.loss_hparams, stage, self.turn_on_sp)
            #if epoch % 100 == 0:
            #    self._print_grad_status(stage)
        else:
            print(f'WARNING: reached epoch_num ({self.epoch_num}) before {stage} stopped early')
            self.completed = False
        # 結果 (utils.post_process) は最後の段階で検証損失が最良だった重みから作る。保存される .pth と同じ重み
        self._restore_best(stage, epoch)
        self._qz_cache = None # 学習後の推論 (utils.post_process) は x から計算するので GPU 上のキャッシュは解放
        self.epoch_num = epoch
        self.last_stage = stage
        print(f'end {stage} at epoch {epoch}')

    def _schedule_done(self, epoch, stage):
        # その段階の温度が 1 まで下がり、段階 1 では warm-up も終わったか。それまでは段階を終える判定と学習率の引き下げをしない
        # (温度が下がる途中は学習の目的関数が変わっていき、warm-up 中は学習率が小さい。最良の重みの記録はこの間も行う)
        temp = self.loss_hparams[{'train_x': 'kld_temp', 'train_pi': 'binom_temp', 'train_sp': 'sp_temp'}[stage]]
        done = temp.temp <= temp.last_temp
        if stage == 'train_x':
            done = done and epoch >= self.warmup_epochs
        return done

    def _restore_best(self, stage, epoch):
        # 段階の終わりに、その段階で検証損失が最良だった重み (StageChange が保存) に戻す。
        # その段階で一度も評価しなかった場合は今の重みを保存し、保存されるモデルと結果を揃える。
        # 各段階の (名前, 最良のエポック, 終えたエポック) を stage_log に残す (workflow が used_parameters.txt に書く)
        self.stage_log.append((stage, self.stagemanager.best_epoch, epoch))
        if self.stagemanager.saved:
            self.model.load_state_dict(torch.load(self.save_model_path, map_location=self.device))
        else:
            torch.save(self.model.state_dict(), self.save_model_path)



    def get_forward(self, s):
        z, d, qz, qd, px_z_ld, pu_zd_ld = self.model(s)
        return(z, d, qz, qd, px_z_ld, pu_zd_ld)

    def init_optimizer(self, lr):
        self.optimizer = torch.optim.Adam(self.model.parameters(), lr=lr)
    
    def update_optimizer_lr(self, lr):
        for param_group in self.optimizer.param_groups:
            param_group['lr'] = lr
        
        # スケジューラを更新する場合(patience/factor は __init__ で受けた値を引き継ぐ)
        if self.lr_scheduler != None:
            self.lr_scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(self.optimizer, 'min', factor=self.lr_factor, patience=self.lr_patience)

    def _sample_and_save_z(self, epoch, stage):
            """
            全データセットの潜在変数zをサンプリングして保存する
            """
            print(f"\nエポック {epoch}: 全データのzをサンプリングして保存しています...")
            self.model.eval()  # 評価モードに設定
            
            # 保存先のサブディレクトリを作成
            save_dir = os.path.dirname(self.save_model_path)
            z_trace_dir = os.path.join(save_dir, 'z_trace')
            os.makedirs(z_trace_dir, exist_ok=True)
            
            epoch_dir = os.path.join(z_trace_dir, f'epoch_{epoch}_{stage}')
            os.makedirs(epoch_dir, exist_ok=True)
            
            # メタデータを記録
            metadata = {
                'epoch': epoch,
                'stage': stage,
                'dataset_sizes': {},
                'z_dim': self.model.enc_z.h2mu.out_features
            }
            
            with torch.no_grad():  # 勾配計算を無効化
                # 全データセットに対してzをサンプリング
                # edm.x, edm.xnorm_mat, edm.d, edm.bには全データが含まれています
                
                # 1. トレーニングデータ
                train_indices = self.edm.train_idx
                train_size = len(train_indices)
                process_batch_size = min(300, train_size)  # メモリ効率のためのバッチサイズ
                
                all_train_z = []
                all_train_z_mean = []
                all_sample_info = []
                
                # バッチ処理
                for start_idx in range(0, train_size, process_batch_size):
                    end_idx = min(start_idx + process_batch_size, train_size)
                    batch_indices = train_indices[start_idx:end_idx]
                    
                    # バッチデータの取得
                    x = self.edm.x[batch_indices].to(self.device)
                    x_norm_mat = self.edm.xnorm_mat[batch_indices].to(self.device)
                    d = self.edm.d[batch_indices].to(self.device)
                    b = self.edm.b[batch_indices].to(self.device)
                    sample_info = self.edm.sample_info[batch_indices].to(self.device)
                    
                    # モデルの前方伝播
                    z, qz, qz_mu, xld, theta_x, phi, pi_z, phi_0, phi_k, alpha = self.model(x, sample_info)
                    
                    # CPU上に移動してNumPy配列に変換
                    all_train_z.append(z.cpu().numpy())
                    all_train_z_mean.append(qz_mu.cpu().numpy())
                    all_sample_info.append(sample_info.cpu().numpy())
                    
                    # メモリ解放
                    del z, qz, qz_mu, xld, theta_x, phi, pi_z, phi_0, phi_k, alpha
                    torch.cuda.empty_cache()
                    
                    if (start_idx + process_batch_size) % 1000 == 0 or end_idx == train_size:
                        print(f"  トレーニングデータ: {end_idx}/{train_size} サンプル処理済み")
                
                # NumPy配列を連結して保存
                train_z = np.concatenate(all_train_z, axis=0) if all_train_z else np.array([])
                train_z_mean = np.concatenate(all_train_z_mean, axis=0) if all_train_z_mean else np.array([])
                
                np.save(os.path.join(epoch_dir, 'train_z.npy'), train_z)
                np.save(os.path.join(epoch_dir, 'train_z_mean.npy'), train_z_mean)
                metadata['dataset_sizes']['train'] = len(train_z)
                
                # 2. 検証データとテストデータ
                for dataset_name, indices in [('val', self.edm.val_idx), ('test', self.edm.test_idx)]:
                    dataset_size = len(indices)
                    process_batch_size = min(300, dataset_size)
                    
                    all_z = []
                    all_z_mean = []
                    
                    for start_idx in range(0, dataset_size, process_batch_size):
                        end_idx = min(start_idx + process_batch_size, dataset_size)
                        batch_indices = indices[start_idx:end_idx]
                        
                        # バッチデータの取得
                        x = self.edm.x[batch_indices].to(self.device)
                        x_norm_mat = self.edm.xnorm_mat[batch_indices].to(self.device)
                        d = self.edm.d[batch_indices].to(self.device)
                        b = self.edm.b[batch_indices].to(self.device)
                        sample_info = self.edm.sample_info[batch_indices].to(self.device)
                        
                        # モデルの前方伝播
                        z, qz, qz_mu, xld, theta_x, phi, pi_z, phi_0, phi_k, alpha = self.model(x, sample_info)
                        
                        # CPU上に移動してNumPy配列に変換
                        all_z.append(z.cpu().numpy())
                        all_z_mean.append(qz_mu.cpu().numpy())
                        
                        # メモリ解放
                        del z, qz, qz_mu, xld, theta_x, phi, pi_z, phi_0, phi_k, alpha
                        torch.cuda.empty_cache()
                    
                    # NumPy配列を連結して保存
                    dataset_z = np.concatenate(all_z, axis=0) if all_z else np.array([])
                    dataset_z_mean = np.concatenate(all_z_mean, axis=0) if all_z_mean else np.array([])
                    
                    np.save(os.path.join(epoch_dir, f'{dataset_name}_z.npy'), dataset_z)
                    np.save(os.path.join(epoch_dir, f'{dataset_name}_z_mean.npy'), dataset_z_mean)
                    metadata['dataset_sizes'][dataset_name] = len(dataset_z)
                    
                    print(f"  {dataset_name}データ: {dataset_size}/{dataset_size} サンプル処理済み")
                    
                # 3. 全データセット
                all_z = []
                all_z_mean = []
                all_sample_info = []
                total_size = len(self.edm.x)
                process_batch_size = min(300, total_size)
                
                for start_idx in range(0, total_size, process_batch_size):
                    end_idx = min(start_idx + process_batch_size, total_size)
                    
                    # バッチデータの取得
                    x = self.edm.x[start_idx:end_idx].to(self.device)
                    x_norm_mat = self.edm.xnorm_mat[start_idx:end_idx].to(self.device)
                    d = self.edm.d[start_idx:end_idx].to(self.device)
                    b = self.edm.b[start_idx:end_idx].to(self.device)
                    sample_info = self.edm.sample_info[start_idx:end_idx].to(self.device)
                    
                    # モデルの前方伝播
                    z, qz, qz_mu, xld, theta_x, phi, pi_z, phi_0, phi_k, alpha = self.model(x, sample_info)
                    
                    # CPU上に移動してNumPy配列に変換
                    all_z.append(z.cpu().numpy())
                    all_z_mean.append(qz_mu.cpu().numpy())
                    all_sample_info.append(sample_info.cpu().numpy())
                    
                    # メモリ解放
                    del z, qz, qz_mu, xld, theta_x, phi, pi_z, phi_0, phi_k, alpha, sample_info
                    torch.cuda.empty_cache()
                    
                    if (start_idx + process_batch_size) % 1000 == 0 or end_idx == total_size:
                        print(f"  全データ: {end_idx}/{total_size} サンプル処理済み")
                
                # NumPy配列を連結して保存
                all_data_z = np.concatenate(all_z, axis=0)
                all_data_z_mean = np.concatenate(all_z_mean, axis=0)
                all_data_sample_info = np.concatenate(all_sample_info, axis=0) if all_sample_info else np.array([])
                
                np.save(os.path.join(epoch_dir, 'all_z.npy'), all_data_z)
                np.save(os.path.join(epoch_dir, 'all_z_mean.npy'), all_data_z_mean)
                np.save(os.path.join(epoch_dir, 'sample_info.npy'), all_data_sample_info)
                metadata['dataset_sizes']['all'] = len(all_data_z)

                def convert_numpy_types(obj):
                    if isinstance(obj, np.integer) or isinstance(obj, np.int64):
                        return int(obj)
                    elif isinstance(obj, np.floating) or isinstance(obj, np.float32) or isinstance(obj, np.float64):
                        return float(obj)
                    elif isinstance(obj, np.ndarray):
                        return obj.tolist()
                    elif isinstance(obj, dict):
                        return {k: convert_numpy_types(v) for k, v in obj.items()}
                    elif isinstance(obj, list):
                        return [convert_numpy_types(item) for item in obj]
                    else:
                        return obj
                    
                # メタデータをJSON形式で保存
                import json
                converted_metadata = convert_numpy_types(metadata)
                with open(os.path.join(epoch_dir, 'metadata.json'), 'w') as f:
                    json.dump(converted_metadata, f, indent=2)
                
                # エポックリストを追跡
                if not hasattr(self, 'z_trace_epochs'):
                    self.z_trace_epochs = []
                self.z_trace_epochs.append(epoch)
                
                # エポックリストをファイルに書き込み
                with open(os.path.join(z_trace_dir, 'z_trace_epochs.txt'), 'w') as f:
                    for e in self.z_trace_epochs:
                        f.write(f"{e},{stage}\n")
                
                print(f"エポック {epoch} のzを保存しました: {epoch_dir}")

    def evaluate_z_clustering(self, epoch, stage, resolutions=[0.2, 0.5, 1.0, 2.0]):
        """
        サンプリングした潜在変数zに対してLeidenクラスタリングを実行し、評価指標を計算する
        
        Parameters:
        -----------
        epoch : int
            現在のエポック数
        stage : str
            現在の学習段階 ('train_x', 'train_pi', 'train_sp')
        resolutions : list of float
            Leidenクラスタリングの解像度パラメータのリスト
        """
        import scanpy as sc
        import anndata
        import os
        import igraph as ig
        from sklearn.metrics import silhouette_score
        import matplotlib.pyplot as plt
        import numpy as np
        from scipy import sparse
        
        print(f"\nエポック {epoch}: zのクラスタリング評価を実行中...")
        
        # 保存先のディレクトリパス
        save_dir = os.path.dirname(self.save_model_path)
        z_trace_dir = os.path.join(save_dir, 'z_trace')
        epoch_dir = os.path.join(z_trace_dir, f'epoch_{epoch}_{stage}')
        cluster_dir = os.path.join(epoch_dir, 'clustering')
        os.makedirs(cluster_dir, exist_ok=True)
        
        # zデータを読み込む
        z_all = np.load(os.path.join(epoch_dir, 'all_z_mean.npy'))
        print('z_all.shape:', z_all.shape)
        sample_info = np.load(os.path.join(epoch_dir, 'sample_info.npy'))
        
        # AnnDataオブジェクトを作成
        print('sample_info.shape:', sample_info.shape)
        adata_z = anndata.AnnData(X=z_all, obs={'sample_info': sample_info[:,0]})
        
        # 近傍グラフの構築
        sc.pp.neighbors(adata_z, use_rep='X')
        
        # Leidenクラスタリングの実行（複数の解像度で）
        results = {}
        
        for res in resolutions:
            print(f"  解像度 {res} でクラスタリング中...")
            # クラスタリング実行
            sc.tl.leiden(adata_z, resolution=res, key_added=f'leiden_res{res}')
            
            # モジュラリティの計算（igraphを使用）
            import igraph as ig
            from sklearn.neighbors import kneighbors_graph
            
            # Get connectivity graph from AnnData object
            if 'neighbors' in adata_z.uns and 'connectivities' in adata_z.uns['neighbors']:
                g = adata_z.uns['neighbors']['connectivities']
            elif hasattr(adata_z, 'obsp') and 'connectivities' in adata_z.obsp:
                g = adata_z.obsp['connectivities']
            else:
                print("  Warning: 'connectivities' not found, creating basic adjacency graph")
                g = kneighbors_graph(adata_z.X, n_neighbors=10, mode='distance')
                g = 0.5 * (g + g.T)  # Make symmetric
            sources, targets = g.nonzero()
            weights = g[sources, targets]
            if isinstance(weights, np.matrix):
                weights = weights.A1
            
            # igraphグラフを構築
            g_ig = ig.Graph(n=g.shape[0], edges=list(zip(sources, targets)), directed=False, edge_attrs={'weight': weights})
            
            # クラスタIDを数値に変換
            clusters = adata_z.obs[f'leiden_res{res}'].astype('category').cat.codes.values
            
            # モジュラリティを計算
            modularity = g_ig.modularity(clusters, weights='weight')
            
            # シルエット係数の計算
            if len(np.unique(adata_z.obs[f'leiden_res{res}'])) > 1:  # クラスタが2つ以上ある場合のみ
                try:
                    # カテゴリ型のクラスタラベルを数値に変換
                    cluster_labels = adata_z.obs[f'leiden_res{res}'].astype('category').cat.codes.values
                    silhouette = silhouette_score(z_all, cluster_labels)
                except Exception as e:
                    print(f"    シルエット係数の計算中にエラーが発生: {e}")
                    silhouette = float('nan')
            else:
                silhouette = float('nan')
            
            # 結果を保存
            results[str(res)] = {
                'modularity': float(modularity),
                'silhouette': float(silhouette),
                'n_clusters': int(len(np.unique(adata_z.obs[f'leiden_res{res}'])))
            }
        
        # 結果をJSONとして保存
        def convert_numpy_types(obj):
            if isinstance(obj, np.integer) or isinstance(obj, np.int64):
                return int(obj)
            elif isinstance(obj, np.floating) or isinstance(obj, np.float32) or isinstance(obj, np.float64):
                return float(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            elif isinstance(obj, dict):
                return {k: convert_numpy_types(v) for k, v in obj.items()}
            elif isinstance(obj, list):
                return [convert_numpy_types(item) for item in obj]
            else:
                return obj
                
        import json
        converted_results = convert_numpy_types(results)
        with open(os.path.join(cluster_dir, 'clustering_metrics.json'), 'w') as f:
            json.dump(converted_results, f, indent=2)
        
        # TensorBoardにも記録（利用可能な場合）
        if self.tb_writer is not None:
            for res, metrics in results.items():
                self.tb_writer.add_scalar(f'clustering/modularity_res{res}', metrics['modularity'], epoch)
                self.tb_writer.add_scalar(f'clustering/silhouette_res{res}', metrics['silhouette'], epoch)
                self.tb_writer.add_scalar(f'clustering/n_clusters_res{res}', metrics['n_clusters'], epoch)
        
        # UMAPも計算して可視化
        try:
            print("  UMAPを計算中...")
            sc.tl.umap(adata_z)
            
            # 各解像度のクラスタリング結果をプロット
            for res in resolutions:
                print(f"  解像度 {res} のクラスタリング結果をプロット中...")
                plt.figure(figsize=(10, 8))
                sc.pl.umap(adata_z, color=f'leiden_res{res}', show=False)
                plt.savefig(os.path.join(cluster_dir, f'umap_leiden_res{res}.png'), dpi=300)
                plt.close()
                
            plt.figure(figsize=(10, 8))
            sc.pl.umap(adata_z, color='sample_info', show=False)
            plt.savefig(os.path.join(cluster_dir, 'umap_sample_info.png'), dpi=300)
            plt.close()
        except Exception as e:
            print(f"  UMAPの計算中にエラーが発生しました: {e}")
        
        # AnnDataオブジェクトを保存
        adata_z.write(os.path.join(cluster_dir, 'adata_z_clustered.h5ad'))
        
        print(f"エポック {epoch} の z クラスタリング評価を完了しました")
        return results
