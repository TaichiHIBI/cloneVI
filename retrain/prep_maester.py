"""MAESTER preprocessing for the model input: raw counts of the 2,000 most variable genes, and the 26 variants analysed
in the original study (Miller et al. 2022)."""
import os

import scanpy as sc

ROOT = os.environ.get('CLONEVI_ROOT', os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')))
MUTATIONS = ['2593', '683', '6205', '9164', '3628', '6185', '9321', '15140', '779', '10158', '6293', '15812', '11864', '1323',
             '1171', '1222', '3988', '1415', '9753', '6243', '8697', '15299', '6218', '2623', '4104', '2241']


def load():
    adata = sc.read_h5ad(f'{ROOT}/data/h5ad_files_RNA/rna_exp.h5ad')
    adata = adata.raw.to_adata()
    adata.layers['raw'] = adata.X
    adata_mg = sc.read_h5ad(f'{ROOT}/data/h5ad_files_RNA/rna_mgv_counts.h5ad')
    adata_mg = adata_mg[:, MUTATIONS]

    adata.obs.rename(columns={"cluster_id": "genotype", "seurat_clusters": "lsi_clusters"}, inplace=True)
    adata_mg.obs.rename(columns={"cluster_id": "genotype"}, inplace=True)
    adata.obs.reset_index(inplace=True, names=['barcode'])
    adata_mg = adata_mg[adata_mg.obs['barcode'].isin(adata.obs['barcode']), :]
    adata = adata[adata.obs['barcode'].isin(adata_mg.obs['barcode']), :]
    adata.obs.reset_index(inplace=True, drop=True)
    adata = adata[adata.obs.sort_values(by='barcode').index].copy()
    adata.obs.reset_index(inplace=True, drop=True)
    adata_mg.obs.reset_index(inplace=True, drop=True)

    adata.obs['genotype'] = adata.obs['MT_clone_num'].astype('int').astype('category')
    adata.obs['MT_clone'] = adata.obs['MT_clone'].astype('category')
    adata_mg.obs['genotype'] = adata.obs['genotype']
    adata.var_names_make_unique()
    adata_mg.var_names_make_unique()
    adata.obs['CellType'] = adata.obs['CellType'].astype('category')

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata = adata[:, adata.var.highly_variable]
    adata.X = adata.layers['raw']
    return adata, adata_mg
