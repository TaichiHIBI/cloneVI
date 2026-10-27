#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N lung_04pk
#$ -l cp -l s_vmem=4G -l mem_req=4G
#$ -o logs/04_pk.$JOB_ID.$TASK_ID.out
#$ -e logs/04_pk.$JOB_ID.$TASK_ID.err
set -euo pipefail
# アレイ: 最小セットを tar -> results_<GSM>.tar.gz。完了後に中間物を掃除(DATA_POLICY)。
WORKDIR=$HOME/GSE302113
SHEET=${SHEET:-$WORKDIR/samples_lung.tsv}
cd "$WORKDIR"

read SRR SAMPLE < <(awk -v n="$SGE_TASK_ID" "NR==n" "$SHEET")
[ -n "${SAMPLE:-}" ] || { echo "ERROR: sheet行取得失敗 task=$SGE_TASK_ID" >&2; exit 1; }
echo "[task $SGE_TASK_ID] SAMPLE=$SAMPLE"

TAR="$WORKDIR/results_${SAMPLE}.tar.gz"
if [ -f "$TAR" ]; then echo "[skip] 既に tar あり: $TAR"; exit 0; fi
OUT="$WORKDIR/$SAMPLE/outs"
EXPORT="$WORKDIR/mgatk_out/$SAMPLE/export"
[ -d "$OUT" ] || { echo "ERROR: cellranger outs無し $OUT" >&2; exit 1; }
[ -f "$EXPORT/coverage.mtx" ] || { echo "ERROR: mgatk export無し $EXPORT" >&2; exit 1; }

PACK="$WORKDIR/pack_$SAMPLE"; rm -rf "$PACK"; mkdir -p "$PACK/cellranger" "$PACK/mgatk"
cp -r "$OUT/filtered_peak_bc_matrix" "$PACK/cellranger/"
for x in peaks.bed summary.csv singlecell.csv fragments.tsv.gz fragments.tsv.gz.tbi web_summary.html; do
  [ -e "$OUT/$x" ] && cp "$OUT/$x" "$PACK/cellranger/" || true
done
cp -r "$EXPORT/"* "$PACK/mgatk/"

tar czf "$TAR" -C "$WORKDIR" "pack_$SAMPLE"
echo "DONE [$SAMPLE] -> $TAR ($(du -sh "$TAR" | cut -f1))"

# --- クリーンアップ(DATA_POLICY) : tar済の中間物を削除 ---
rm -rf "$PACK"
rm -rf "$WORKDIR/$SAMPLE"                 # cellranger outs(bamは03で削除済)
rm -rf "$WORKDIR/mgatk_out/$SAMPLE"       # mgatk final/temp/export
echo "cleanup done [$SAMPLE]; 残置: $(basename "$TAR")"
