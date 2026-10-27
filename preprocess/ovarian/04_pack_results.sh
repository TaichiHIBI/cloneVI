#!/bin/bash
# このPCに戻す最小セットを tar 化（巨大な FASTQ / possorted_bam.bam は残す）
set -euo pipefail
WORKDIR=$HOME/GSE302113
SAMPLE=GSM9096478
OUT="$WORKDIR/$SAMPLE/outs"
PACK="$WORKDIR/pack_$SAMPLE"
rm -rf "$PACK"; mkdir -p "$PACK/cellranger" "$PACK/mgatk"

# --- cellranger: Signac/h5ad に必要な最小限 ---
cp -r "$OUT/filtered_peak_bc_matrix" "$PACK/cellranger/"
cp "$OUT/peaks.bed"                  "$PACK/cellranger/" 2>/dev/null || true
cp "$OUT/summary.csv"                "$PACK/cellranger/" 2>/dev/null || true
cp "$OUT/singlecell.csv"             "$PACK/cellranger/" 2>/dev/null || true
# Signac は fragments とその index を要求（やや大きいが ~1GB 程度）
cp "$OUT/fragments.tsv.gz"           "$PACK/cellranger/"
cp "$OUT/fragments.tsv.gz.tbi"       "$PACK/cellranger/" 2>/dev/null || true

# --- mgatk: 行列 export 一式 ---
cp -r "$WORKDIR/mgatk_out/export/"*  "$PACK/mgatk/"

cd "$WORKDIR"
tar czf "results_${SAMPLE}.tar.gz" -C "$WORKDIR" "pack_$SAMPLE"
echo "DONE -> $WORKDIR/results_${SAMPLE}.tar.gz"
du -sh "results_${SAMPLE}.tar.gz"
echo "これを GSE302113/prototype_GSM9096478/ に転送してください。"
