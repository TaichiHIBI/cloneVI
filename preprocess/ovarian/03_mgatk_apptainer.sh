#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N gse_mgatk
#$ -pe def_slot 32
#$ -l s_vmem=16G -l mem_req=16G
#$ -o logs/03_mgatk.$JOB_ID.out -e logs/03_mgatk.$JOB_ID.err
set -euo pipefail
# ============================================================================
# mtscATAC mtDNA variant calling : cellranger-atac BAM -> mgatk tenx -> MTX
# 旧 03_mgatk_apptainer.sh(mgatk_v1.sif=v0.5.3 + bcall, 1コア5h) からの改訂版。
#   - biocontainer mgatk v0.7.0 を使用 (mgatk + R/SummarizedExperiment + picard同梱)
#   - bcall(細胞毎にBAM分割=実質1コア) -> tenx(コア数チャンク並列) に変更
#   - 罠対策: --cleanenv + PYTHONNOUSERSITE=1 でホスト ~/.local の古いmgatk漏れ込み防止
#   - 罠対策: MALLOC_ARENA_MAX/BLASスレッド=1 で仮想メモリ膨張を抑制
#             (実RSSは数G。shirokaneは s_vmem==mem_req 必須なので 16G/slot 確保)
#   コマンドは ~/mgatk.sh / 03_mgatk.sh の tenx 用法を踏襲。
# ============================================================================
NSLOTS=${NSLOTS:-32}

WORKDIR=$HOME/GSE302113
SAMPLE=GSM9096478
APPT=$(ls -d /usr/local/package/apptainer/*/bin/apptainer 2>/dev/null | head -1)
SIF=$WORKDIR/mgatk_bioc_0.7.0.sif
MITO=hg38                                  # cellranger-atac GRCh38 = chrM

cd "$WORKDIR"; mkdir -p logs
[ -x "$APPT" ] || { echo "ERROR: apptainer無し" >&2; exit 1; }
[ -f "$SIF" ]  || { echo "ERROR: SIF無し: $SIF" >&2; exit 1; }

OUT="$WORKDIR/$SAMPLE/outs"
BAM="$OUT/possorted_bam.bam"
[ -f "$BAM" ] || { echo "ERROR: BAM無し: $BAM" >&2; exit 1; }

zcat -f "$OUT"/filtered_peak_bc_matrix/barcodes.tsv* > "$WORKDIR/filtered_barcodes.txt"
echo "[info] barcodes: $(wc -l < "$WORKDIR/filtered_barcodes.txt") 個 / cores: $NSLOTS"

BIND="-B /home -B /hshare1"; [ -d /swork ] && BIND="$BIND -B /swork"
ENVS="--env PYTHONNOUSERSITE=1 --env PYTHONPATH= --env MALLOC_ARENA_MAX=2 \
 --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 --env NUMEXPR_NUM_THREADS=1"
RUN() { "$APPT" exec --cleanenv $ENVS $BIND "$SIF" "$@"; }

echo "[sanity] chrM idxstats:"
RUN python -c "import pysam; print(\"\n\".join(l for l in pysam.idxstats(\"$BAM\").splitlines() if l.split()[0] in (\"chrM\",\"MT\")))" || true

echo "[run] mgatk tenx ..."
RUN mgatk tenx -i "$BAM" -n "$SAMPLE" -o "$WORKDIR/mgatk_out" \
  -g "$MITO" -bt CB -b "$WORKDIR/filtered_barcodes.txt" -c "$NSLOTS"

RDS=$(ls "$WORKDIR/mgatk_out/final/"*.rds 2>/dev/null | grep -v signac | head -1 || true)
[ -n "$RDS" ] || RDS=$(ls "$WORKDIR/mgatk_out/final/"*.rds 2>/dev/null | head -1 || true)
[ -n "$RDS" ] || { echo "ERROR: RDS無し。mgatk出力を確認" >&2; exit 1; }
echo "[run] RDS=$RDS -> mtx export"

EXPORT_R="$WORKDIR/03b_export_mgatk_mtx.R"
[ -f "$EXPORT_R" ] || EXPORT_R="$WORKDIR/shirokane/03b_export_mgatk_mtx.R"
[ -f "$EXPORT_R" ] || { echo "ERROR: export R無し" >&2; exit 1; }
RUN Rscript "$EXPORT_R" "$RDS" "$WORKDIR/mgatk_out/export"

echo "DONE."; ls -la "$WORKDIR/mgatk_out/export"
