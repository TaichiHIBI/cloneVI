import torch
import numpy as np
import utils
import random
import pandas as pd
from scipy import sparse as sp


# ---- メモリ最適化版 (2026-06) -------------------------------------------------
# 大規模 ATAC (cells x peaks ~ 1e5 x 1e5) を全 dense 化すると count_mat / norm_mat /
# 一時コピーで 100GB 超 → OOM。対策として:
#   * count(x) は scipy 疎行列のまま保持し、__getitem__ でミニバッチ分だけ dense 化
#   * norm_mat は rank-1 (scale * row_mean[i] * col_mean[j]) なので行列を作らず
#     バッチ時に再構成（元の値と数学的に同一）
# exp.py / model 側は dense バッチを受け取る形のまま無変更で動く。
# ----------------------------------------------------------------------------

def _as_csr(x):
    return x.tocsr() if sp.issparse(x) else sp.csr_matrix(x)


class VaeDataSet(torch.utils.data.Dataset):
    def __init__(self, x, row_mean, col_mean, scale, d, b, sample_info=None, transform=None, pre_transform=None):
        self.x = x                    # scipy csr (cells x peaks) 疎のまま保持
        self.row_mean = row_mean      # torch (cells,)
        self.col_mean = col_mean      # torch (peaks,)
        self.scale = float(scale)
        self.d = d
        self.b = b
        self.sample_info = sample_info

    def __len__(self):
        return(self.x.shape[0])

    def __getitem__(self, idx):
        idx_x = torch.from_numpy(np.asarray(self.x[idx].todense()).ravel().astype(np.float32))
        # xnorm は rank-1 (scale*row_mean[i]*col_mean[j])。dense 展開(85859次元)して転送するのは
        # col_mean を毎細胞ぶん重複転送する無駄なので、row_mean[i] スカラーだけ返し GPU で再構成する。
        idx_row_mean = self.row_mean[idx]
        idx_d = self.d[idx]
        idx_b = self.b[idx]
        idx_sample_info = self.sample_info[idx]
        return(idx_x, idx_row_mean, idx_d, idx_b, idx_sample_info)


def worker_init_fn(worker_id):
    np.random.seed(42 + worker_id)
    random.seed(42 + worker_id)
    torch.manual_seed(42 + worker_id)
    torch.cuda.manual_seed(42 + worker_id)
    torch.cuda.manual_seed_all(42 + worker_id)

class VaeDataManager():
    def __init__(self, adata, adata_mg, batch_size=300, num_workers=2, val_frac=0.05, test_frac=0.1, seed=42, sample_info=None):
        self.split_seed = seed  # train/val/test 分割の seed (モデルの seed とは別。workflow が記録する)
        # ---- count(x): 疎のまま保持し rank-1 norm の係数だけ算出 ----
        X = _as_csr(adata.X)                     # cells x peaks, 疎
        N, M = X.shape
        row_mean = torch.from_numpy(np.asarray(X.sum(axis=1)).ravel().astype(np.float32) / M)  # (cells,)
        col_mean = torch.from_numpy(np.asarray(X.sum(axis=0)).ravel().astype(np.float32) / N)  # (peaks,)
        overall_mean = float(X.sum()) / (N * M)
        # norm_mat[i,j] = overall_mean * (row_mean[i]*col_mean[j]) / mean(row_mean*col_mean)
        #   mean over all (i,j) of row_mean[i]*col_mean[j] = mean(row_mean)*mean(col_mean)
        denom = float(row_mean.mean()) * float(col_mean.mean())
        scale = overall_mean / denom if denom != 0 else 0.0

        # ---- mtSNV: d, b は変異数が小さい(conf5~1.5k)ので dense で保持 ----
        d = torch.Tensor(adata_mg.X.toarray()) if type(adata_mg.X) != np.ndarray else torch.Tensor(adata_mg.X)
        b = torch.Tensor(adata_mg.layers['all_mut'].toarray()) if type(adata_mg.layers['all_mut']) != np.ndarray else torch.Tensor(adata_mg.layers['all_mut'])
        print('adata_mg_shape', adata_mg.X.shape)
        if sample_info is not None:
            sample_info = adata.obs[sample_info]
            print("Converting DataFrame to tensor")
            sample_info = torch.Tensor(pd.get_dummies(sample_info, dtype=int).values)
            print('sample_info_shape:', sample_info.shape)
        else:
            sample_info = torch.zeros(N, 1)

        # 参照用(疎/係数)。フル dense の self.x/self.xnorm_mat は持たない(メモリ削減)
        self.x = X
        self.row_mean = row_mean
        self.col_mean = col_mean
        self.scale = scale
        self.d = d
        self.b = b
        self.sample_info = sample_info

        # Create full dataset (疎保持、__getitem__ でバッチ密化)
        full_dataset = VaeDataSet(X, row_mean, col_mean, scale, d, b, sample_info)

        # Calculate dataset sizes
        total_num = N
        self.test_num = int(total_num * test_frac)
        self.val_num = int(total_num * val_frac)
        self.train_num = total_num - self.val_num - self.test_num
        print(total_num)

        # Use random_split instead of RandomState
        train_dataset, val_dataset, test_dataset = torch.utils.data.random_split(
            full_dataset,
            [self.train_num, self.val_num, self.test_num],
            generator=torch.Generator().manual_seed(seed)
        )

        # Extract indices for compatibility with existing code
        self.train_idx = train_dataset.indices
        self.val_idx = val_dataset.indices
        self.test_idx = test_dataset.indices

        # val/test は小さい部分集合なので一度だけ dense 化（norm は rank-1 で再構成）
        def _dense_x(rows):
            return torch.from_numpy(np.asarray(X[rows].todense()).astype(np.float32))
        def _dense_norm(rows):
            rm = row_mean[rows].reshape(-1, 1)
            return (scale * rm) * col_mean.reshape(1, -1)

        # Extract validation data
        self.val_x = _dense_x(self.val_idx)
        self.val_xnorm_mat = _dense_norm(self.val_idx)
        self.val_d = d[self.val_idx]
        self.val_b = b[self.val_idx]
        self.val_sample_info = sample_info[self.val_idx] if sample_info is not None else None

        # Extract test data
        self.test_x = _dense_x(self.test_idx)
        self.test_xnorm_mat = _dense_norm(self.test_idx)
        self.test_d = d[self.test_idx]
        self.test_b = b[self.test_idx]
        self.test_sample_info = sample_info[self.test_idx] if sample_info is not None else None

        # Use train_dataset directly for DataLoader
        self.dataset = train_dataset
        self.train_loader = torch.utils.data.DataLoader(
            self.dataset, batch_size=batch_size, num_workers=num_workers, shuffle=True, worker_init_fn=worker_init_fn)
