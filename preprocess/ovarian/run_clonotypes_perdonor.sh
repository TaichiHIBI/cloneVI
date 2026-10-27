#!/usr/bin/env bash
# ドナー単位 mtSNV 検出 + FindClonotypes を **ドナーごとに別Rプロセス** で逐次実行し
# （メモリを完全解放しながら）、最後に per-donor 出力をコホート単位ファイルに結合する。
# 使い方: bash run_clonotypes_perdonor.sh <COHORT>     COHORT=Ovarian|Lung
set -euo pipefail
COHORT="${1:-Ovarian}"
CELLS="cohort_${COHORT}_atac_cells.csv"
ENVNAME="Signac"

donors=$(tail -n +2 "$CELLS" | cut -d, -f3 | tr -d '"' | sort -u)
echo "COHORT=$COHORT donors:"; echo "$donors" | sed 's/^/  /'

for d in $donors; do
  vs="cohort_${COHORT}_${d}_variant_stats.tsv"
  cl="cohort_${COHORT}_${d}_clonotype.csv"
  if [[ -s "$vs" && -s "$cl" ]]; then
    echo "[skip] $d (出力あり)"
    continue
  fi
  echo "===== $d 実行 ($(date +%H:%M:%S)) ====="
  conda run -n "$ENVNAME" Rscript clonotypes_perdonor_signac.R "$COHORT" "$d"
  echo "----- $d 完了 ($(date +%H:%M:%S)) -----"
done

# --- 結合: per-donor → cohort ---
MERGED_VS="cohort_${COHORT}_perdonor_variant_stats.tsv"
MERGED_CL="cohort_${COHORT}_clonotypes.csv"
first=1
: > "$MERGED_VS"
for d in $donors; do
  vs="cohort_${COHORT}_${d}_variant_stats.tsv"
  [[ -s "$vs" ]] || { echo "WARN: $vs 無し（$d はスキップ?）"; continue; }
  if [[ $first -eq 1 ]]; then cat "$vs" > "$MERGED_VS"; first=0
  else tail -n +2 "$vs" >> "$MERGED_VS"; fi
done
first=1
: > "$MERGED_CL"
for d in $donors; do
  cl="cohort_${COHORT}_${d}_clonotype.csv"
  [[ -s "$cl" ]] || continue
  if [[ $first -eq 1 ]]; then cat "$cl" > "$MERGED_CL"; first=0
  else tail -n +2 "$cl" >> "$MERGED_CL"; fi
done
echo "MERGED -> $MERGED_VS ($(wc -l < "$MERGED_VS") rows), $MERGED_CL ($(wc -l < "$MERGED_CL") rows)"
echo "high.conf union variants: $(tail -n +2 "$MERGED_VS" | awk -F'\t' '$NF=="TRUE"{print $3}' | sort -u | wc -l)"
