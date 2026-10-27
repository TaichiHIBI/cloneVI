#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N gse_sra
#$ -pe def_slot 8
#$ -l s_vmem=8G
#$ -l mem_req=8G
#$ -o logs/01_sra.$JOB_ID.out
#$ -e logs/01_sra.$JOB_ID.err
set -euo pipefail
NSLOTS=${NSLOTS:-8}   # qsub時はUGEが設定。bashで直接実行する時の既定。
# 総mem = s_vmem(8G) x def_slot(8) = 64G。prefetchは1コア/網I/O、fasterq-dumpは$NSLOTS並列。
# 速度の効きどころは fasterq-dump の並列と temp 領域の速さ。slotは8で十分(増やしても頭打ち)。

module use /usr/local/package/modulefiles/
module load sra-toolkit
module load pigz

# ===== CONFIG (要編集) =====
WORKDIR=$HOME/GSE302113
SRR=SRR34441260
SAMPLE=GSM9096478
# aria2c の置き場所(.bashrc不要・ここでPATHに通す)。違う場所なら変更。
ARIA2_DIR=$HOME/opt/aria2
[ -d "$ARIA2_DIR" ] && export PATH="$ARIA2_DIR:$PATH"
# fasterq-dump は 2–3× の一時容量(数百GB)が必要。NFS home だと激遅 → 高速・大容量scratchを推奨。
# 例: SCRATCH=/swork/<group>/<user>/gse_tmp  （ノードローカル $TMPDIR があればそれが最速だが容量に注意）
SCRATCH=${TMPDIR:-$WORKDIR/tmp}
# ===========================
mkdir -p "$WORKDIR/logs" "$SCRATCH" "$WORKDIR/fastq_raw"
cd "$WORKDIR"

# [1] DL: aria2c(並列, 速い) があれば S3公開ミラーから、無ければ prefetch にフォールバック。
SRA_FILE="sra/$SRR/$SRR.sra"
if command -v aria2c >/dev/null 2>&1; then
  echo "[1] aria2c 並列DL (S3 mirror, 16接続)"
  # 静的aria2cはCA束を同梱しないので、システムのCAを指定。無ければ検証スキップ(公開データ)。
  CA=""
  for c in /etc/pki/tls/certs/ca-bundle.crt /etc/ssl/certs/ca-certificates.crt /etc/ssl/certs/ca-bundle.crt; do
    [ -f "$c" ] && { CA="$c"; break; }
  done
  if [ -n "$CA" ]; then TLS=(--ca-certificate="$CA"); else TLS=(--check-certificate=false); fi
  echo "    TLS: ${TLS[*]}"
  aria2c -x16 -s16 -k50M -c "${TLS[@]}" -o "$SRR.sra" \
    "https://sra-pub-run-odp.s3.amazonaws.com/sra/$SRR/$SRR"
  SRA_FILE="$WORKDIR/$SRR.sra"
else
  echo "[1] prefetch $SRR  (aria2c が PATH に無いため。~/bin に静的aria2cを入れると速い)"
  prefetch --max-size 200G -O sra "$SRR"
fi

echo "[2] fasterq-dump (10x: テクニカル読み含む, $NSLOTS スレッド, temp=$SCRATCH)"
fasterq-dump --split-files --include-technical -e "$NSLOTS" -t "$SCRATCH" \
  -O fastq_raw "$SRA_FILE"

echo "[3] read長を判定して 10x(R1=genomic, R2=16bp barcode, R3=genomic) にマッピング"
mkdir -p fastq
declare -A LEN
for f in fastq_raw/${SRR}_*.fastq; do
  L=$(awk 'NR%4==2{s+=length($0); n++} n>=2000{exit} END{print int(s/n)}' "$f")
  LEN[$f]=$L
  echo "    $f  modal_len=$L"
done
# 16bp(±2) を barcode(R2)、長い方2つを R1,R3（元番号順）に割当
R2=""; GEN=()
for f in $(ls fastq_raw/${SRR}_*.fastq | sort); do
  L=${LEN[$f]}
  if [ "$L" -ge 14 ] && [ "$L" -le 18 ]; then R2="$f"; else GEN+=("$f"); fi
done
if [ -z "$R2" ] || [ "${#GEN[@]}" -lt 2 ]; then
  echo "!! 自動判定に失敗。LEN一覧を見て手動で R1/R2(16bp)/R3 を割当ててください。" >&2
  exit 1
fi
R1="${GEN[0]}"; R3="${GEN[1]}"
echo "    -> R1=$R1  R2(barcode)=$R2  R3=$R3"

ln -sf "$WORKDIR/$R1" "fastq/${SAMPLE}_S1_L001_R1_001.fastq"
ln -sf "$WORKDIR/$R2" "fastq/${SAMPLE}_S1_L001_R2_001.fastq"
ln -sf "$WORKDIR/$R3" "fastq/${SAMPLE}_S1_L001_R3_001.fastq"
echo "[4] gzip (並列)"
if command -v pigz >/dev/null 2>&1; then
  # pigz: 各ファイルを $NSLOTS スレッドで圧縮（逐次でも十分速い）
  for r in R1 R2 R3; do
    pigz -1 -p "$NSLOTS" -c "fastq/${SAMPLE}_S1_L001_${r}_001.fastq" > "fastq/${SAMPLE}_S1_L001_${r}_001.fastq.gz"
  done
else
  # pigz が無ければ 3 本を同時(バックグラウンド)に gzip → 約3倍
  for r in R1 R2 R3; do
    gzip -1 -c "fastq/${SAMPLE}_S1_L001_${r}_001.fastq" > "fastq/${SAMPLE}_S1_L001_${r}_001.fastq.gz" &
  done
  wait
fi
rm -f fastq/${SAMPLE}_S1_L001_R1_001.fastq fastq/${SAMPLE}_S1_L001_R2_001.fastq fastq/${SAMPLE}_S1_L001_R3_001.fastq
echo "DONE. cellranger-atac 用 FASTQ: $WORKDIR/fastq/  (--sample=$SAMPLE)"
echo "※ R2 が 16bp barcode であることを必ず確認してから 02 を実行。"
