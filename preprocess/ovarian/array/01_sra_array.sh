#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N lung_01sra
#$ -pe def_slot 8
#$ -l s_vmem=8G -l mem_req=8G
#$ -o logs/01_sra.$JOB_ID.$TASK_ID.out
#$ -e logs/01_sra.$JOB_ID.$TASK_ID.err
set -euo pipefail
# アレイ: samples_lung.tsv の SGE_TASK_ID 行を処理。
# 生fastq(213G級)はノードローカル $TMPDIR に出し、gz だけ Lustre fastq/<GSM>/ へ。
# reads=0 のSRRは fasterq-dump が空 -> 自動的に exit 1(タスク失敗)。
NSLOTS=${NSLOTS:-8}
WORKDIR=$HOME/GSE302113
SHEET=${SHEET:-$WORKDIR/samples_lung.tsv}
cd "$WORKDIR"

read SRR SAMPLE < <(awk -v n="$SGE_TASK_ID" "NR==n" "$SHEET")
[ -n "${SRR:-}" ] && [ -n "${SAMPLE:-}" ] || { echo "ERROR: sheet行取得失敗 task=$SGE_TASK_ID" >&2; exit 1; }
echo "[task $SGE_TASK_ID] SRR=$SRR SAMPLE=$SAMPLE host=$(hostname)"

OUTDIR="$WORKDIR/fastq/$SAMPLE"
if ls "$OUTDIR"/${SAMPLE}_S1_L001_R3_001.fastq.gz >/dev/null 2>&1; then
  echo "[skip] 既に fastq あり: $OUTDIR"; exit 0
fi
mkdir -p "$OUTDIR" logs

module use /usr/local/package/modulefiles/
module load sra-toolkit
module load pigz 2>/dev/null || echo "[info] pigz無し -> gzip使用"
command -v fasterq-dump >/dev/null 2>&1 || { echo "ERROR: fasterq-dump無し(sra-toolkit module確認)" >&2; exit 1; }

LOCAL="${TMPDIR:-/tmp}/sra_$SAMPLE"; mkdir -p "$LOCAL/raw"
ARIA2_DIR=$HOME/opt/aria2; [ -d "$ARIA2_DIR" ] && export PATH="$ARIA2_DIR:$PATH"

# [1] DL -> ノードローカル
SRA="$LOCAL/$SRR.sra"
if command -v aria2c >/dev/null 2>&1; then
  echo "[1] aria2c DL (S3 mirror)"
  CA=""; for c in /etc/pki/tls/certs/ca-bundle.crt /etc/ssl/certs/ca-certificates.crt; do [ -f "$c" ] && { CA="$c"; break; }; done
  if [ -n "$CA" ]; then TLS=(--ca-certificate="$CA"); else TLS=(--check-certificate=false); fi
  aria2c -x16 -s16 -k50M -c "${TLS[@]}" -d "$LOCAL" -o "$SRR.sra" \
    "https://sra-pub-run-odp.s3.amazonaws.com/sra/$SRR/$SRR"
else
  echo "[1] prefetch"; prefetch --max-size 200G -O "$LOCAL" "$SRR"; SRA="$LOCAL/$SRR/$SRR.sra"
fi

# [2] fasterq-dump -> ノードローカル
echo "[2] fasterq-dump (temp/out=$LOCAL)"
fasterq-dump --split-files --include-technical -e "$NSLOTS" -t "$LOCAL" -O "$LOCAL/raw" "$SRA"

# [3] read長で R1/R2(16bp)/R3 を割当
declare -A LEN; for f in "$LOCAL"/raw/${SRR}_*.fastq; do
  [ -e "$f" ] || { echo "ERROR: fastq無し(reads=0の可能性) SRR=$SRR" >&2; exit 1; }
  LEN[$f]=$(awk "NR%4==2{s+=length(\$0);n++} n>=2000{exit} END{print int(s/n)}" "$f"); echo "    $f len=${LEN[$f]}"
done
R2=""; GEN=()
for f in $(ls "$LOCAL"/raw/${SRR}_*.fastq | sort); do
  L=${LEN[$f]}; if [ "$L" -ge 14 ] && [ "$L" -le 18 ]; then R2="$f"; else GEN+=("$f"); fi
done
[ -n "$R2" ] && [ "${#GEN[@]}" -ge 2 ] || { echo "ERROR: R1/R2/R3 自動判定失敗 SRR=$SRR" >&2; exit 1; }
R1="${GEN[0]}"; R3="${GEN[1]}"; echo "    R1=$R1 R2=$R2 R3=$R3"

# [4] 圧縮 -> Lustre fastq/<GSM>/  (pigz優先, 無ければgzip 3並列)
echo "[4] compress -> $OUTDIR"
declare -A SRCF=( [R1]="$R1" [R2]="$R2" [R3]="$R3" )
if command -v pigz >/dev/null 2>&1; then
  for r in R1 R2 R3; do pigz -1 -p "$NSLOTS" -c "${SRCF[$r]}" > "$OUTDIR/${SAMPLE}_S1_L001_${r}_001.fastq.gz"; done
else
  for r in R1 R2 R3; do gzip -1 -c "${SRCF[$r]}" > "$OUTDIR/${SAMPLE}_S1_L001_${r}_001.fastq.gz" & done
  wait
fi
# 完成検証(空でないこと)
for r in R1 R2 R3; do
  [ -s "$OUTDIR/${SAMPLE}_S1_L001_${r}_001.fastq.gz" ] || { echo "ERROR: $r 圧縮失敗/空" >&2; exit 1; }
done
rm -rf "$LOCAL"
echo "DONE [$SAMPLE]: $(ls -la "$OUTDIR" | grep -c fastq.gz) gz"
