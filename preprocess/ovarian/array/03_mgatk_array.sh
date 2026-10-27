#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N lung_03mg
#$ -pe def_slot 16
#$ -l s_vmem=16G -l mem_req=16G
#$ -o logs/03_mg.$JOB_ID.$TASK_ID.out
#$ -e logs/03_mg.$JOB_ID.$TASK_ID.err
set -euo pipefail
# アレイ: mgatk tenx (biocontainer)。確定レシピ: cleanenv+PYTHONNOUSERSITE, arena/BLAS抑制, s_vmem=mem_req=16G。
# 出力は mgatk_out/<GSM>/。完了後 possorted_bam.bam を削除(DATA_POLICY)。
NSLOTS=${NSLOTS:-16}
WORKDIR=$HOME/GSE302113
SHEET=${SHEET:-$WORKDIR/samples_lung.tsv}
APPT=$(ls -d /usr/local/package/apptainer/*/bin/apptainer 2>/dev/null | head -1)
SIF=$WORKDIR/mgatk_bioc_0.7.0.sif
MITO=hg38
cd "$WORKDIR"

read SRR SAMPLE < <(awk -v n="$SGE_TASK_ID" "NR==n" "$SHEET")
[ -n "${SAMPLE:-}" ] || { echo "ERROR: sheet行取得失敗 task=$SGE_TASK_ID" >&2; exit 1; }
echo "[task $SGE_TASK_ID] SAMPLE=$SAMPLE host=$(hostname) cores=$NSLOTS"

OUT="$WORKDIR/mgatk_out/$SAMPLE"
if [ -f "$OUT/export/coverage.mtx" ]; then echo "[skip] 既にmgatk export あり: $OUT"; exit 0; fi
rm -rf "$OUT"   # requeue対策: 前回中断の残骸を一掃(checkGrep誤判定を防ぐ)
BAM="$WORKDIR/$SAMPLE/outs/possorted_bam.bam"
[ -f "$BAM" ] || { echo "ERROR: BAM無し $BAM" >&2; exit 1; }
[ -x "$APPT" ] && [ -f "$SIF" ] || { echo "ERROR: apptainer/SIF無し" >&2; exit 1; }

BC="$WORKDIR/$SAMPLE/outs/filtered_peak_bc_matrix/barcodes.tsv"
BCLIST="$WORKDIR/mgatk_out/$SAMPLE.barcodes.txt"; mkdir -p "$WORKDIR/mgatk_out"
zcat -f "$BC"* > "$BCLIST"
echo "[info] barcodes: $(wc -l < "$BCLIST")"

BIND="-B /home -B /hshare1"; [ -d /swork ] && BIND="$BIND -B /swork"
ENVS="--env PYTHONNOUSERSITE=1 --env PYTHONPATH= --env MALLOC_ARENA_MAX=2 \
 --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 --env NUMEXPR_NUM_THREADS=1"
RUN() { "$APPT" exec --cleanenv $ENVS $BIND "$SIF" "$@"; }

RUNCWD="$WORKDIR/mgatk_out/.snk_$SAMPLE"; rm -rf "$RUNCWD"; mkdir -p "$RUNCWD"; cd "$RUNCWD"   # snakemakeロックをサンプル毎に分離
echo "[run] mgatk tenx -> $OUT"
RUN mgatk tenx -i "$BAM" -n "$SAMPLE" -o "$OUT" -g "$MITO" -bt CB -b "$BCLIST" -c "$NSLOTS"

RDS=$(ls "$OUT/final/"*.rds 2>/dev/null | grep -v signac | head -1 || true)
[ -n "$RDS" ] || RDS=$(ls "$OUT/final/"*.rds 2>/dev/null | head -1 || true)
[ -n "$RDS" ] || { echo "ERROR: RDS無し" >&2; exit 1; }
RUN Rscript "$WORKDIR/03b_export_mgatk_mtx.R" "$RDS" "$OUT/export"
[ -f "$OUT/export/coverage.mtx" ] || { echo "ERROR: export失敗" >&2; exit 1; }

rm -f "$BCLIST"
# bam は消費済 -> 削除(DATA_POLICY)
rm -f "$BAM" "$BAM.bai"
cd "$WORKDIR"; rm -rf "$RUNCWD"
echo "DONE [$SAMPLE]: export $(ls "$OUT/export"/*.mtx | wc -l) mtx; bam削除"
