library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(Matrix)
library(SummarizedExperiment)

rna_exp <- readRDS("MAESTER/clonalhemato/BPDCN712_Seurat_with_TCR.rds")
rna_mgv <- readRDS("MAESTER/clonalhemato/TenX_BPDCN712_All_mr3_maegatk.rds")

if(!dir.exists('MAESTER/rds_output/')){
  dir.create('MAESTER/rds_output/')
}

for (n in 1:length(assayNames(rna_mgv))) {
  writeMM(assay(rna_mgv,n), paste0("MAESTER/rds_output/", assayNames(rna_mgv)[n], ".mtx"))
}

write(colnames(rna_mgv),"MAESTER/rds_output/barcodes.tsv")
write.table(rowRanges(rna_mgv),"MAESTER/rds_output/genes.tsv",sep='\t')

if(!dir.exists('MAESTER/h5ad_files_RNA/')){
  dir.create('MAESTER/h5ad_files_RNA/')
}

# クローン情報を番号に変換。NAを0に変換
rna_exp$MT_clone_num <- is.na(rna_exp$MT_clone_num)
rna_exp$MT_clone_num[is.na(rna_exp$MT_clone_num)] <- 0
rna_exp$MT_clone <- as.character(rna_exp$MT_clone)

SaveH5Seurat(rna_exp, 'MAESTER/h5ad_files_RNA/rna_exp')
Convert('MAESTER/h5ad_files_RNA/rna_exp.h5seurat', dest='h5ad')
