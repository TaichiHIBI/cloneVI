import scanpy as sc
import torch
import dataset, module, utils, exp, plot
import os
import datetime
from torch.utils.tensorboard import SummaryWriter
import numpy as np
import pandas as pd

def resolve_model_params(adata, model_params, mutations, sample_info=None):
    # 既定値を補い、入力データの大きさから次元を決める (estimate_probability と qz_from_checkpoint で共通)。model_params をその場で書き換えて返す
    default_params = {
        'x_dim': None,
        'xz_dim': 10,
        'enc_z_h_dim': 50,
        'dec_z_h_dim': 50,
        'num_enc_z_layers': 2,
        'num_dec_z_layers': 2,
        'k': 20,
        'use_phi_common': False,  # バックグラウンド変異がある場合はTrueに
        'freeze_z_after_x': True,  # False で train_pi/sp でも enc_z を解凍（z へ勾配を流す）
        'cache_qz_after_x': True  # enc_z 凍結時、train_pi/sp は q(z|x) をキャッシュして x を読まない (False で旧来の経路)
    }
    # Initialize model_params if None
    if model_params is None:
        model_params = {}

    # Merge default parameters with provided parameters
    for key, value in default_params.items():
        if key not in model_params:
            model_params[key] = value

    # Set dimensions based on input data
    model_params['x_dim'] = adata.shape[1]
    model_params['xz_dim'] = np.sqrt(np.sqrt(adata.shape[1]).astype(int)).astype(int) if np.sqrt(np.sqrt(adata.shape[1]).astype(int)) > model_params['xz_dim'] else model_params['xz_dim']
    model_params['enc_z_h_dim'] = np.sqrt(adata.shape[1]).astype(int) if np.sqrt(adata.shape[1]).astype(int) > model_params['enc_z_h_dim'] else model_params['enc_z_h_dim']
    model_params['dec_z_h_dim'] = np.sqrt(adata.shape[1]).astype(int) if np.sqrt(adata.shape[1]).astype(int) > model_params['dec_z_h_dim'] else model_params['dec_z_h_dim']
    model_params['mutations'] = mutations
    #model_params['n_samples'] = sample_info.shape[1] if sample_info is not None else 1
    model_params['n_samples'] = adata.obs[sample_info].nunique() if sample_info is not None else 1
    return(model_params)


def qz_from_checkpoint(adata, model_path, mutations, model_params=None, sample_info=None, batch_size=1024):
    # 保存済みの学習の .pth から、全細胞 (adata の行順) の q(z|x) の平均と標準偏差を計算する。
    # adata・mutations・model_params・sample_info はその学習と同じものを渡す (同じ形のモデルを作って重みを読む)。
    # 段階 2 以降は encoder を凍結しているので、どの段階で保存された .pth でも encoder は段階 1 の終わりのもの。
    # 戻り値を estimate_probability(..., qz=...) に渡すと、段階 1 を学習し直さずに段階 2 から始められる。
    # qz を与えて学習した run の .pth は encoder を学習していないので使えない (その run の qz.npz を使う)
    mp = resolve_model_params(adata, dict(model_params or {}), mutations, sample_info)
    for key in ('freeze_z_after_x', 'cache_qz_after_x'): # VaeExperiment が使う設定で scVAE の引数ではない
        mp.pop(key, None)
    model = module.scVAE(**mp)
    model.load_state_dict(torch.load(model_path, map_location='cpu'))
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model.to(device).eval()
    X = adata.X
    N = X.shape[0]
    # sample_info の one-hot は VaeDataManager と同じ作り方
    si = torch.Tensor(pd.get_dummies(adata.obs[sample_info], dtype=int).values) if sample_info is not None else torch.zeros(N, 1)
    mus, scales = [], []
    with torch.no_grad():
        for s0 in range(0, N, batch_size):
            s1 = min(s0 + batch_size, N)
            xb = X[s0:s1]
            xb = np.asarray(xb.todense()) if hasattr(xb, 'todense') else np.asarray(xb)
            mu, scale = model.encode_qz(torch.from_numpy(xb.astype(np.float32)).to(device), si[s0:s1].to(device))
            mus.append(mu.cpu().numpy()); scales.append(scale.cpu().numpy())
    return(np.concatenate(mus), np.concatenate(scales))


def estimate_probability(
        adata, adata_mg, mutations, save_dir, loss_hparams,
        model_params = None,
        epoch_num = 20000,
        lr = 4e-3, lr_scheduler = True, lr_patience = 10, lr_factor = 0.1,
        stage_patience_x = 10, stage_patience_pi = 20,
        val_ratio = 0.05, test_ratio = 0.1, batch_size = 300, num_workers = 2,
        n_neighbors=30, min_dist=0.5, celltype='celltype', genotype='genotype', train_mode = 'train_all',
        sample_info=None, use_tensorboard = False, experiment_name=None, warmup_epochs=3000, max_trainx_epochs=None, seed=42, reconst_mat=False,
        qz=None, qz_source=None):
    # qz: (平均, 標準偏差) の組 (各 細胞数 x 潜在次元、adata の行順)。与えると段階 1 を学習せず、この q(z|x) から段階 2 を始める。
    #     qz_from_checkpoint の出力 (保存済みの段階 1 を使い回す) や、外部の埋め込み (標準偏差 None なら z を平均に固定) を渡す。
    # qz_source: qz の出どころ (used_parameters.txt に記録する文字列)

    utils.set_seed(seed)
    model_params = resolve_model_params(adata, model_params, mutations, sample_info)
    if qz is not None:
        qz_mu, qz_scale = qz
        model_params['xz_dim'] = np.asarray(qz_mu).shape[1]
        train_mode = 'train_pi'

    # save model path
    dt_now = datetime.datetime.now()
    now_month=dt_now.strftime('%Y%m')
    now_day=dt_now.strftime('%d')
    now_time_detail=dt_now.strftime('%H%M%S')
    if experiment_name == None:
        experiment_name = input('Enter experiment name: ')
    save_path = os.path.join(save_dir, experiment_name, now_month, now_day, now_time_detail)
    save_model_dir = os.path.join(save_path, 'models')
    save_model_path = os.path.join(save_model_dir, 'bestmodel.pth')
    save_plt_path = os.path.join(save_path, 'plots')
    os.makedirs(save_path, exist_ok=True)
    os.makedirs(save_plt_path, exist_ok=True)
    os.makedirs(save_model_dir, exist_ok=True)

    if use_tensorboard == True:
        tb_writer = SummaryWriter(log_dir=os.path.join(save_path, 'logs'))
    else:
        tb_writer = None

    # Save parameters to a file
    params_log = {
        'adata_shape': adata.shape,
        'adata_mg_shape': adata_mg.shape,
        'model_params': model_params,
        'loss_hparams': {k: v.temp if 'temp' in k else v for k, v in loss_hparams.items()},
        'epoch_num': epoch_num,
        'lr': lr,
        'val_ratio': val_ratio,
        'test_ratio': test_ratio,
        'batch_size': batch_size,
        'num_workers': num_workers,
        'n_neighbors': n_neighbors,
        'min_dist': min_dist,
        'celltype_label': celltype,
        'genotype_label': genotype,
        'train_mode': train_mode,
        'use_multisamples': sample_info is not None,
        'use_sparse_term': loss_hparams.get('lmd_pi', 0) != 0 or loss_hparams.get('lmd_phi', 0) != 0,
        'use_tensorboard': use_tensorboard,
        'experiment_name': experiment_name,
        'warmup_epochs': warmup_epochs,
        'max_trainx_epochs': max_trainx_epochs,
        'stage_patience_x': stage_patience_x,
        'stage_patience_pi': stage_patience_pi,
        'lr_patience': lr_patience,
        'lr_factor': lr_factor,
        'seed': seed,
        'qz_source': qz_source if qz is not None else None
    }
    
    # Save parameters to a text file
    with open(os.path.join(save_path, 'used_parameters.txt'), 'w') as f:
        for key, value in params_log.items():
            f.write(f'{key}: {value}\n')
            print(f'{key}: {value}')
    # load model
    
    cloneVI_exp = exp.VaeExperiment(adata, adata_mg, model_params, lr, lr_scheduler, val_ratio, test_ratio, batch_size, num_workers, loss_hparams=loss_hparams, save_model_path=save_model_path, epoch_num=epoch_num, sample_info=sample_info, tb_writer=tb_writer, warmup_epochs=warmup_epochs, max_trainx_epochs=max_trainx_epochs,
                                    stage_patience_x=stage_patience_x, stage_patience_pi=stage_patience_pi, lr_patience=lr_patience, lr_factor=lr_factor)
    # 分割を記録する (分割の seed は上の seed とは別で VaeDataManager の既定値)。段階 2 から学習し直すときにも使う
    edm = cloneVI_exp.edm
    with open(os.path.join(save_path, 'used_parameters.txt'), 'a') as f:
        f.write(f'split_seed: {edm.split_seed}\n')
    np.savez(os.path.join(save_path, 'split_idx.npz'), train=np.asarray(edm.train_idx), val=np.asarray(edm.val_idx), test=np.asarray(edm.test_idx))
    if qz is not None:
        cloneVI_exp.set_qz_cache(qz_mu, qz_scale)
        # この run の .pth の encoder は学習していない (初期値のまま) ので、結果は .pth ではなく与えた q(z|x) から再現する。
        # qz_from_checkpoint にこの run の .pth を渡さず、ここに保存した平均・標準偏差を使う
        np.savez(os.path.join(save_path, 'qz.npz'), mu=np.asarray(qz_mu), **({} if qz_scale is None else {'scale': np.asarray(qz_scale)}))

    cloneVI_exp.train_total(train_mode)

    if use_tensorboard == True:
        cloneVI_exp.tb_writer.close()

    # train_total は最後の段階で検証損失が最良だった重みに戻して終わるので、結果の h5ad と models/ の .pth は同じ重み。
    # epoch_num に達して打ち切られた場合などは completed: False
    with open(os.path.join(save_path, 'used_parameters.txt'), 'a') as f:
        f.write(f'completed: {cloneVI_exp.completed}\nlast_stage: {cloneVI_exp.last_stage}\nlast_epoch: {cloneVI_exp.epoch_num}\n')
        f.write(f'stages (stage, best epoch, end epoch): {cloneVI_exp.stage_log}\n')

    adata_mg, test_adata_mg = utils.post_process(adata, adata_mg, cloneVI_exp, n_neighbors, min_dist,
                                                 qz_mu=None if qz is None else qz_mu)

    pl_list = ['loss_plots', 'embedding_plots', 'pi_plots', 'pi_learn_plots', 'phi_plots', 'emb_pi_plots']
    for i in pl_list:
        os.makedirs(os.path.join(save_plt_path, i), exist_ok=True)
    #plot.π_learn_plot(cloneVI_exp.epoch_num, cloneVI_exp.model.pi_list, save_plt_path, cloneVI_exp.test_loss_list[-1])
    plot.π_plot(adata_mg.obsm['pi_sorted'], save_plt_path)
    plot.embedding_plot(adata_mg.obsm['X_cloneVI_umap'], adata, adata_mg, save_plt_path, celltype, genotype)
    plot.pi_emb_plot(adata_mg.obsm['X_cloneVI_umap'], adata_mg.obsm['pi_sorted'], save_plt_path)
    plot.phi_plot(adata_mg.varm['phi_sorted'], adata_mg.varm['phi_0'], adata_mg.varm['phi_k'], adata_mg.uns['alpha'], save_plt_path, mutations)

    adata_mg.write(os.path.join(save_path, 'cloneVI_result.h5ad'))
    test_adata_mg.write(os.path.join(save_path, 'cloneVI_test_result.h5ad'))
    # change model name
    save_rename_path = os.path.join(save_path, f'models/{cloneVI_exp.test_loss_list[-1]}.pth')
    #os.remove(save_model_path)
    os.rename(save_model_path, save_rename_path)
    return(adata, adata_mg)