import os
import sys
import pandas as pd
import numpy as np
import umap
import seaborn as sns
from matplotlib import pyplot as plt
import scipy
import scipy.stats
import scipy.sparse
import torch.distributions as dist
import torch
import math
import random
from torch.optim.lr_scheduler import _LRScheduler

def set_seed(seed=42):
    """再現性のためにシードを設定する"""
    random.seed(seed)  # Pythonの組み込み乱数のシード
    np.random.seed(seed)  # NumPyのシード
    torch.manual_seed(seed)  # PyTorchのシード
    torch.cuda.manual_seed(seed)  # CUDA のシード（GPU使用時）
    torch.cuda.manual_seed_all(seed)  # 複数GPU使用時のシード
    torch.backends.cudnn.deterministic = True  # CuDNNの決定論的動作を有効化
    torch.backends.cudnn.benchmark = False  # CuDNNの最適化を無効化（完全再現のため）


class StageChange:
    """earlystoppingクラス"""

    def __init__(self, path, patience=5):
        """引数：最小値の非更新数カウンタ、表示設定、モデル格納path"""

        self.patience = patience    #設定ストップカウンタ
        self.counter = 0            #現在のカウンタ値
        self.best_score = None      #ベストスコア
        self.flag = False   #ストップフラグ
        self.val_loss_min = np.inf  #前回のベストスコア記録用
        self.path = path
        self.saved = False  # この段階の最良の重みを path に保存したか (段階の終わりに読み戻す。exp.VaeExperiment._restore_best)
        self.best_epoch = None

    @staticmethod
    def criterion(loss_list, stage):
        # 各段階の検証損失。温度で重みを付けない値なので、温度が下がる途中の評価とも比べられる
        # (段階 1 は旧来 loss_list[0] = KL/温度 + x を使っていた。2026-10-10 に KL + x に変更)
        if stage == 'train_x':
            return loss_list[3] + loss_list[4]
        elif stage == 'train_pi':
            return loss_list[5]
        elif stage == 'train_sp':
            return loss_list[5] + loss_list[6] + loss_list[7]
        print('invalid stage name.')
        sys.exit(1)

    def __call__(self, loss_list, model, stage, epoch=None, can_stop=True):
        """
        特殊(call)メソッド
        実際に学習ループ内で最小lossを更新したか否かを計算させる部分。
        最良の重みは段階のすべての評価から選んで保存する。止める判定 (flag) は can_stop のときだけ
        (warm-up と温度の下降が終わるまでは止めない。exp.VaeExperiment._schedule_done)
        """
        val_loss = self.criterion(loss_list, stage)
        score = -val_loss
        if self.best_score is None or score >= self.best_score:  # 最初の評価、またはベストスコアを更新した場合
            self.best_score = score
            self.best_epoch = epoch
            self.checkpoint(val_loss, model)
            self.counter = 0
        else:  # ベストスコアを更新できなかった場合
            self.counter += 1
            if can_stop and self.counter >= self.patience:
                self.flag = True

    def checkpoint(self, val_loss, model):
        torch.save(model.state_dict(), self.path) # ベストモデルをpathに保存
        self.saved = True

# Temperature scheduler
class DecreasingTemp:
    def __init__(self, dec_temp_steps, init_temp, last_temp=1):
        self.temp = init_temp
        self.last_temp = last_temp
        self.dec_temp_rate = math.pow(last_temp / init_temp, 1.0 / dec_temp_steps)

    def step(self):
        self.temp = max(self.dec_temp_rate * self.temp, self.last_temp)

    @property
    def value(self):
        return self.temp
    
class DecreasingTemp_adapt:
    def __init__(self, dec_temp_steps, init_temp, last_temp=1, monitor_ratio=False):
        self.temp = init_temp
        self.init_temp = init_temp
        self.last_temp = last_temp
        self.dec_temp_rate = math.pow(last_temp / init_temp, 1.0 / dec_temp_steps)
        self.monitor_ratio = monitor_ratio
        self.prev_loss = None
        
    def step(self, current_loss=None):
        if self.monitor_ratio and current_loss is not None:
            # 損失が増加した場合、温度を少し戻す
            if self.prev_loss is not None and current_loss > self.prev_loss:
                self.temp = min(self.temp / self.dec_temp_rate, self.init_temp)
            else:
                self.temp = max(self.temp * self.dec_temp_rate, self.last_temp)
            self.prev_loss = current_loss
        else:
            # 通常の温度減少
            self.temp = max(self.dec_temp_rate * self.temp, self.last_temp)
    @property
    def value(self):
        return self.temp
    
def temp_scheduler(params, stage, turn_on_sp, current_loss=False):
    if stage == 'train_x':
        params['kld_temp'].step()
    elif turn_on_sp == True:
        if stage == 'train_pi':
            params['binom_temp'].step()
        elif stage == 'train_sp' and current_loss != False:
            params['sp_temp'].step(current_loss)
        else:
            params['sp_temp'].step()
    else:
        if stage == 'train_pi':
            params['binom_temp'].step()
        elif stage == 'train_sp':
            params['sp_temp'] = 0
        else:
            print('invalid stage name.')
            sys.exit(1)
        

# sorting and merge
def sorting_func(pi, phi, phi0, phik, alpha, k):
    re_clone = {(i): pi.mean(axis=0).cpu()[i] for i in range(k)}
    re_clone_sorted = sorted(re_clone.items(), key=lambda x:x[1], reverse=True)
    re_clone_sorted = [re_clone_sorted[i][0] for i in range(k)]
    pi_sort = torch.Tensor(pi.size(dim=0), k)
    for i in range(k):
      pi_sort[:,i] = pi[:,re_clone_sorted[i]]

    phi_sort = torch.Tensor(k, phi.size(dim=1))
    phi0_sort = torch.Tensor(k, phi.size(dim=1))
    phik_sort = torch.Tensor(k, phi.size(dim=1))
    alpha_sort = torch.Tensor(k, phi.size(dim=1))

    for i in range(k):
        phi_sort[i] = phi[re_clone_sorted[i]]
        phi0_sort[i] = phi0[re_clone_sorted[i]]
        phik_sort[i] = phik[re_clone_sorted[i]]
        alpha_sort[i] = alpha[re_clone_sorted[i]]

    return(pi_sort, phi_sort, phi0_sort, phik_sort, alpha_sort)

def cal_differ(phi_sort):
    diff = torch.zeros(phi_sort.size(dim=0), phi_sort.size(dim=0))
    for i in range(phi_sort.size(dim=0)):
        for j in range(phi_sort.size(dim=0)):
            tmp = phi_sort[i,:] - phi_sort[j,:]
            tmp = torch.linalg.norm(tmp)
            diff[i,j] = tmp
    return(diff)

# π < 0.01 のものを削除する
def del_clone(pi_sort, phi_sort, diff):
    #index = torch.where((pi_sort.mean(dim=0) > 0.01) | (pi_sort.max(dim=0).values >= 0.2))[0].tolist()
    # Check if pi_sort is a tensor or numpy array
    if isinstance(pi_sort, torch.Tensor):
        index = torch.where((pi_sort.mean(dim=0) > 0.01) | (pi_sort.max(dim=0).values >= 0.2))[0].tolist()
    else:  # numpy array
        index = np.where(np.mean(pi_sort, axis=0) >= 0.01 | (pi_sort.max(dim=0).values >= 0.2))[0].tolist()
    pi_sort = pi_sort[:,index]
    phi_sort = phi_sort[index,:]
    #diff = diff[index[0]:index[-1]+1, index[0]:index[-1]+1] # 少し複雑だが、他と同じことをしている
    return pi_sort, phi_sort, diff

def embed_z(z_mat, n_neighbors=30, min_dist=0.5, seed=42):
    reducer = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, random_state=seed)
    z_embed = reducer.fit_transform(z_mat)
    return(z_embed)

def _batched_infer(model, X, sample_info, batch_size=2048, qz_mu=None):
    """バッチ化＋疎対応の決定的推論（大規模データで OOM 回避）。
    pi は sampled z ではなく qz_mu(事後平均)から計算 → 決定的・バッチ分割で厳密同一。
    大域パラメータ(phi/phi_0/phi_k/alpha)は入力非依存なので model.calc_phi() から取得。
    qz_mu を与えたとき (段階 2 から学習した場合、workflow.estimate_probability の qz) は encoder を通さずその qz_mu を使う。
    戻り値: qz_mu(N,z), pi(N,k), phi, phi_0, phi_k, alpha
    """
    def _dense(Xrows):
        if scipy.sparse.issparse(Xrows):
            return torch.from_numpy(np.asarray(Xrows.todense()).astype(np.float32))
        return torch.as_tensor(Xrows, dtype=torch.float32)
    model.eval()
    N = X.shape[0] if qz_mu is None else len(qz_mu)
    qz_mu_list, pi_list = [], []
    with torch.no_grad():
        phi, phi_0, phi_k, alpha = model.calc_phi()
        for s0 in range(0, N, batch_size):
            s1 = min(s0 + batch_size, N)
            sib = sample_info[s0:s1]
            if qz_mu is None:
                xb = _dense(X[s0:s1])
                mu_b, _ = model.enc_z(torch.cat([xb, sib], axis=1))           # 事後平均(決定的)
            else:
                mu_b = torch.as_tensor(np.asarray(qz_mu[s0:s1]), dtype=torch.float32)
            pi_det = model.dec_z2pi(torch.cat([mu_b, sib], axis=1))        # qz_mu からデコード
            qz_mu_list.append(mu_b)
            pi_list.append(pi_det)
    return torch.cat(qz_mu_list, 0), torch.cat(pi_list, 0), phi, phi_0, phi_k, alpha


def post_process(adata, adata_mg, clonevi_exp, n_neighbors, min_dist, dz_var_prop=0.05, sample_num=10, qz_mu=None):
    # qz_mu: 段階 2 から学習した場合 (workflow.estimate_probability の qz) に与えた全細胞の平均。z と pi はこれから計算する
    clonevi_exp.device = torch.device('cpu')
    clonevi_exp.model = clonevi_exp.model.to(clonevi_exp.device)

    # 全データの推論（バッチ化＋疎対応＋pi は qz_mu から決定的に算出）
    qz_mu_given = qz_mu
    qz_mu, pi, phi, phi_0, phi_k, alpha = _batched_infer(
        clonevi_exp.model, clonevi_exp.edm.x, clonevi_exp.edm.sample_info, qz_mu=qz_mu_given)
    pi_sorted, phi_sorted, phi_0_sorted, phi_k_sorted, alpha_sorted = sorting_func(pi, phi, phi_0, phi_k, alpha, clonevi_exp.model.k)
    pi_deleted, phi_deleted, diff = del_clone(pi_sorted, phi_sorted, cal_differ(phi_sorted))
    z_mat = qz_mu.detach().numpy()
    phi_mat = phi_sorted
    pi_mat = pi_sorted
    phi_0_mat = phi_0_sorted
    phi_k_mat = phi_k_sorted
    alpha_mat = alpha_sorted

    # テストデータの推論（同上）
    test_qz_mu, test_pi, test_phi, test_phi_0, test_phi_k, test_alpha = _batched_infer(
        clonevi_exp.model, clonevi_exp.edm.test_x, clonevi_exp.edm.test_sample_info,
        qz_mu=None if qz_mu_given is None else np.asarray(qz_mu_given)[np.asarray(clonevi_exp.edm.test_idx)])
    test_pi_sorted, test_phi_sorted, test_phi_0_sorted, test_phi_k_sorted, test_alpha_sorted = sorting_func(test_pi, test_phi, test_phi_0, test_phi_k, test_alpha, clonevi_exp.model.k)
    test_pi_deleted, test_phi_deleted, test_diff = del_clone(test_pi_sorted, test_phi_sorted, cal_differ(test_phi_sorted))
    test_z_mat = test_qz_mu.detach().numpy()
    
    # 全データの結果をadata_mgに保存
    adata_mg.obsm['z'] = z_mat
    adata_mg.obsm['pi_sorted'] = pi_mat.detach().numpy()
    adata_mg.varm['phi_sorted'] = phi_mat.T.detach().numpy()
    adata_mg.varm['phi_0'] = phi_0_mat.T.detach().numpy()
    adata_mg.varm['phi_k'] = phi_k_mat.T.detach().numpy()
    adata_mg.uns['alpha'] = alpha_mat.T.detach().numpy()
    adata_mg.obsm['pi_deleted'] = pi_deleted.detach().numpy()
    adata_mg.varm['phi_deleted'] = phi_deleted.T.detach().numpy()
    adata_mg.layers['reconst_VAF'] = (pi_mat @ phi_mat).detach().numpy()

    z_embed = embed_z(z_mat)
    adata_mg.obsm['X_cloneVI_umap'] = z_embed
    
    # テストデータのAnnDataオブジェクトを作成
    test_adata_mg = adata_mg[clonevi_exp.edm.test_idx].copy()
    test_adata_mg.obsm['z'] = test_z_mat
    test_adata_mg.obsm['pi_sorted'] = test_pi_sorted.detach().numpy()
    test_adata_mg.varm['phi_sorted'] = test_phi_sorted.T.detach().numpy()
    test_adata_mg.varm['phi_0'] = test_phi_0_sorted.T.detach().numpy()
    test_adata_mg.varm['phi_k'] = test_phi_k_sorted.T.detach().numpy()
    test_adata_mg.uns['alpha'] = test_alpha_sorted.T.detach().numpy()
    test_adata_mg.obsm['pi_deleted'] = test_pi_deleted.detach().numpy()
    test_adata_mg.varm['phi_deleted'] = test_phi_deleted.T.detach().numpy()
    test_adata_mg.layers['reconst_VAF'] = (test_pi_sorted @ test_phi_sorted).detach().numpy()
    adata_mg.varm['phi_specific'] = test_phi_k_sorted.T.detach().numpy()
    
    test_z_embed = embed_z(test_z_mat)
    test_adata_mg.obsm['X_cloneVI_umap'] = test_z_embed

    return adata_mg, test_adata_mg


def _rank_cols(M):
    """各列を順位に変換する（同点は平均順位）。Spearman = 順位後の Pearson の前段。"""
    try:
        # scipy >= 1.10 は axis 引数をサポート
        return scipy.stats.rankdata(M, axis=0).astype(np.float64)
    except TypeError:
        # 古い scipy 向けフォールバック（pandas の rank は既定で平均順位）
        return pd.DataFrame(M).rank(axis=0, method="average").to_numpy(dtype=np.float64)


def _zscore_cols(R, on_constant="zero"):
    """列ごとに z-score 標準化（母標準偏差 ddof=0）。分散0の列は 0 か NaN にする。"""
    R = R - R.mean(axis=0, keepdims=True)
    s = R.std(axis=0, ddof=0, keepdims=True)
    fill = 0.0 if on_constant == "zero" else np.nan
    return np.divide(R, s, out=np.full_like(R, fill), where=s != 0)


def spearman_cross(A, B, chunk=None, on_constant="zero"):
    """A の各列と B の各列の Spearman 順位相関を総当たりで計算する。

    共有軸は行（例: 細胞）。A は (n, p)、B は (n, q) で列数は異なってよい。
    返り値は (p, q) の相関行列。

    高速化の原理: Spearman は「列ごとに順位変換した後の Pearson」に等しく、
    Pearson は「列ごとに z-score 標準化したベクトルの内積 / n」に等しい。
    よって順位変換は列数分だけで済み、相関本体は単一の行列積に落ちる
    （scipy.stats.spearmanr を組ごとに呼ぶ二重ループより桁違いに速い）。

    Parameters
    ----------
    A, B : array-like, shape (n, p) / (n, q)
        np.ndarray / scipy.sparse / pandas を受ける（内部で密配列に変換）。
        行（共有軸）の数は一致している必要がある。
    chunk : int or None
        B の列をこの本数ずつ分割して処理する。巨大・sparse な B のメモリ抑制用。
    on_constant : {"zero", "nan"}
        分散0の列（全ゼロ遺伝子など）の相関を 0 にするか NaN にするか。
        scipy.stats.spearmanr は NaN を返すので、それに合わせたい場合は "nan"。

    Returns
    -------
    corr : np.ndarray, shape (p, q)
        各 (A列 i, B列 j) の Spearman 順位相関係数。
    """
    A = np.asarray(A.todense()) if scipy.sparse.issparse(A) else np.asarray(A)
    n = A.shape[0]
    if B.shape[0] != n:
        raise ValueError(f"行（共有軸）の数が一致しません: A={n}, B={B.shape[0]}")

    ZA = _zscore_cols(_rank_cols(A), on_constant=on_constant)  # (n, p)

    if chunk is None:
        B_dense = np.asarray(B.todense()) if scipy.sparse.issparse(B) else np.asarray(B)
        ZB = _zscore_cols(_rank_cols(B_dense), on_constant=on_constant)  # (n, q)
        return ZA.T @ ZB / n

    q = B.shape[1]
    out = np.empty((A.shape[1], q), dtype=np.float64)
    for j0 in range(0, q, chunk):
        j1 = min(j0 + chunk, q)
        blk = B[:, j0:j1]
        blk = np.asarray(blk.todense()) if scipy.sparse.issparse(blk) else np.asarray(blk)
        ZB = _zscore_cols(_rank_cols(blk), on_constant=on_constant)
        out[:, j0:j1] = ZA.T @ ZB / n
    return out

