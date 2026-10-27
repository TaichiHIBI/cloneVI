#!/bin/bash
# 肺27本(または引数のシート)を 01->02->03->04 の順次アレイで投入。
# 段間は -hold_jid で全完了バリア。各段 -tc で 128スロット枠内に制御。
#   01: 8core x tc16 =128 / 02,03: 16core x tc8 =128 / 04: cp.q(軽量)
set -euo pipefail
WORKDIR=$HOME/GSE302113
cd "$WORKDIR"; mkdir -p logs
SHEET_REL=${1:-samples_lung.tsv}
SHEET="$WORKDIR/$SHEET_REL"
[ -f "$SHEET" ] || { echo "ERROR: シート無し $SHEET" >&2; exit 1; }
N=$(grep -c . "$SHEET")
echo "SHEET=$SHEET  N=$N 本"
V="SHEET=$SHEET"

j1=$(qsub -terse -t 1-$N -tc 16 -v "$V" array/01_sra_array.sh | cut -d. -f1)
echo "01 SRA-DL    : job $j1  (-t 1-$N -tc16)"
j2=$(qsub -terse -hold_jid "$j1" -t 1-$N -tc 8 -v "$V" array/02_cellranger_array.sh | cut -d. -f1)
echo "02 cellranger: job $j2  (after $j1, -tc8)"
j3=$(qsub -terse -hold_jid "$j2" -t 1-$N -tc 8 -v "$V" array/03_mgatk_array.sh | cut -d. -f1)
echo "03 mgatk     : job $j3  (after $j2, -tc8)"
j4=$(qsub -terse -hold_jid "$j3" -t 1-$N -tc 10 -v "$V" array/04_pack_array.sh | cut -d. -f1)
echo "04 pack      : job $j4  (after $j3)"
echo
echo "進捗: qstat / logs/*.out,err   中断: qdel $j1 $j2 $j3 $j4"
echo "成果物: results_<GSM>.tar.gz (38->ローカル転送対象)"
