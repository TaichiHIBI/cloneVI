#!/bin/bash
# 依存ゼロ(curl + coreutils)の並列セグメントDL。aria2c/conda 不要。
# 高遅延な対米HTTPSを N 本の Range リクエストで並列取得し、結合・サイズ検証する。
# 使い方: bash par_get.sh <URL> <出力ファイル> [並列数=16]
#   例: bash par_get.sh \
#         "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR34441260/SRR34441260" \
#         ~/GSE302113/SRR34441260.sra 16
set -euo pipefail
URL="${1:?usage: par_get.sh URL OUTFILE [NPARTS]}"
OUT="${2:?usage: par_get.sh URL OUTFILE [NPARTS]}"
N="${3:-16}"

SIZE=$(curl -sIL "$URL" | awk 'tolower($1)=="content-length:"{v=$2} END{gsub(/\r/,"",v); print v}')
[ -n "$SIZE" ] || { echo "ERROR: Content-Length 取得失敗 ($URL)" >&2; exit 1; }
curl -sIL "$URL" | grep -qi 'accept-ranges: bytes' || echo "WARN: Range非対応の可能性。続行する。" >&2
echo "size=$SIZE bytes  parts=$N  -> $OUT"

PD="${OUT}.parts"; mkdir -p "$PD"
chunk=$(( (SIZE + N - 1) / N ))
pids=()
for ((i=0; i<N; i++)); do
  s=$(( i*chunk )); e=$(( s+chunk-1 )); (( e>=SIZE )) && e=$(( SIZE-1 ))
  (( s>=SIZE )) && break
  curl -fsSL --retry 5 --retry-delay 3 -r "${s}-${e}" \
       -o "${PD}/part.$(printf '%04d' "$i")" "$URL" &
  pids+=($!)
done
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=1; done
(( fail==0 )) || { echo "ERROR: 一部パートが失敗。再実行してください。" >&2; exit 1; }

cat "${PD}"/part.* > "$OUT"
rm -rf "$PD"
got=$(stat -c %s "$OUT")
echo "downloaded $got / $SIZE bytes"
(( got==SIZE )) || { echo "ERROR: サイズ不一致。$OUT を消して再実行。" >&2; exit 1; }
echo "OK: size一致。先頭16byte:"; head -c 16 "$OUT" | od -c | head -1
echo "次: fasterq-dump --split-files --include-technical -e 8 -t ./tmp -O fastq_raw \"$OUT\""
