#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N lung_02cr
#$ -pe def_slot 16
#$ -l s_vmem=8G -l mem_req=8G
#$ -o logs/02_cr.$JOB_ID.$TASK_ID.out
#$ -e logs/02_cr.$JOB_ID.$TASK_ID.err
set -euo pipefail
# アレイ: cellranger-atac を ノードローカル $TMPDIR 上で実行(Lustre I/O輻輳回避)。
# 完了後、下流に必要な成果物だけ Lustre <GSM>/outs/ にコピー。fastq は消費後削除。
NSLOTS=${NSLOTS:-16}
WORKDIR=$HOME/GSE302113
SHEET=${SHEET:-$WORKDIR/samples_lung.tsv}
CELLRANGER_ATAC=$HOME/opt/cellranger-atac-2.1.0/bin/cellranger-atac
REF=$HOME/opt/refdata-cellranger-arc-GRCh38-2020-A-2.0.0
cd "$WORKDIR"

read SRR SAMPLE < <(awk -v n="$SGE_TASK_ID" "NR==n" "$SHEET")
[ -n "${SAMPLE:-}" ] || { echo "ERROR: sheet行取得失敗 task=$SGE_TASK_ID" >&2; exit 1; }
echo "[task $SGE_TASK_ID] SAMPLE=$SAMPLE host=$(hostname) cores=$NSLOTS"

DEST="$WORKDIR/$SAMPLE/outs"
if [ -f "$DEST/possorted_bam.bam" ]; then echo "[skip] 既に outs あり: $DEST"; exit 0; fi
[ -x "$CELLRANGER_ATAC" ] || { echo "ERROR: cellranger無し" >&2; exit 1; }
FQ="$WORKDIR/fastq/$SAMPLE"
ls "$FQ"/${SAMPLE}_S1_L001_R3_001.fastq.gz >/dev/null 2>&1 || { echo "ERROR: fastq無し $FQ" >&2; exit 1; }

LOCALMEM=$(( NSLOTS * 8 * 95 / 100 ))
RUNROOT="${TMPDIR:-/tmp}/cr_$SAMPLE"; mkdir -p "$RUNROOT"; cd "$RUNROOT"
echo "[run] cellranger-atac on node-local $RUNROOT (localmem=${LOCALMEM}G)"
"$CELLRANGER_ATAC" count \
  --id="$SAMPLE" --reference="$REF" \
  --fastqs="$FQ" --sample="$SAMPLE" \
  --localcores="$NSLOTS" --localmem="$LOCALMEM"

SRC="$RUNROOT/$SAMPLE/outs"
[ -f "$SRC/possorted_bam.bam" ] || { echo "ERROR: cellranger outs生成失敗" >&2; exit 1; }
echo "[copy] 必要物のみ -> $DEST"
mkdir -p "$DEST"
for x in possorted_bam.bam possorted_bam.bam.bai fragments.tsv.gz fragments.tsv.gz.tbi \
         peaks.bed summary.csv singlecell.csv web_summary.html; do
  [ -e "$SRC/$x" ] && cp -a "$SRC/$x" "$DEST/" || echo "  (任意/欠落: $x)"
done
cp -a "$SRC/filtered_peak_bc_matrix" "$DEST/"
rm -rf "$RUNROOT"
# fastq は消費済 -> 削除(DATA_POLICY B)
rm -rf "$FQ"
echo "DONE [$SAMPLE]: $(du -sh "$DEST" | cut -f1) copied; fastq削除"
