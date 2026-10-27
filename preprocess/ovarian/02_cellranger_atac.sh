#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N gse_cratac
#$ -pe def_slot 32
#$ -l s_vmem=8G
#$ -l mem_req=8G
#$ -o logs/02_cratac.$JOB_ID.out
#$ -e logs/02_cratac.$JOB_ID.err
set -euo pipefail
NSLOTS=${NSLOTS:-32}   # qsub時はUGEが設定。bashで直接実行する時の既定。
# 総mem = 8G x 16 = 128G。cellranger-atac は 8コア/64GB が下限・16コア/128GB が推奨。
# さらに速くしたい時は def_slot を 24(→192G) などに上げる(s_vmemは据置)。効果は頭打ち気味。
# --localcores / --localmem は確保量に合わせて自動算出（mem は総量の~95%）。

# cellranger-atac は module に無い(ローカル導入)。実行ファイルとhg38参照をパス指定する。
# ===== CONFIG (要編集) =====
WORKDIR=$HOME/GSE302113
CELLRANGER_ATAC=$HOME/opt/cellranger-atac-2.1.0/bin/cellranger-atac
REF=$HOME/opt/refdata-cellranger-arc-GRCh38-2020-A-2.0.0
SAMPLE=${SAMPLE:-GSM9096478}
RUNID=${RUNID:-$SAMPLE}
# ===========================
mkdir -p "$WORKDIR/logs"
cd "$WORKDIR"
command -v "$CELLRANGER_ATAC" >/dev/null 2>&1 || [ -x "$CELLRANGER_ATAC" ] || { echo "ERROR: cellranger-atac が無い: $CELLRANGER_ATAC" >&2; exit 1; }

LOCALMEM=$(( NSLOTS * 8 * 95 / 100 ))   # GB（総mem=NSLOTS*s_vmem の約95%）
echo "cellranger-atac count ($SAMPLE)  cores=$NSLOTS localmem=${LOCALMEM}G"
"$CELLRANGER_ATAC" count \
  --id="$RUNID" \
  --reference="$REF" \
  --fastqs="$WORKDIR/fastq" \
  --sample="$SAMPLE" \
  --localcores="$NSLOTS" \
  --localmem="$LOCALMEM"

OUT="$WORKDIR/$RUNID/outs"
echo "DONE. 主要出力:"
ls -la "$OUT"/fragments.tsv.gz "$OUT"/peaks.bed "$OUT"/possorted_bam.bam 2>/dev/null || true
ls -la "$OUT"/filtered_peak_bc_matrix/ 2>/dev/null || true
echo "細胞数(filtered barcodes): $(zcat -f "$OUT"/filtered_peak_bc_matrix/barcodes.tsv* | wc -l)"
