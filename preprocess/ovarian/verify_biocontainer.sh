#$ -S /bin/bash
#$ -cwd
#$ -N verify_bioc
#$ -pe def_slot 8
#$ -l s_vmem=8G -l mem_req=8G
#$ -o logs/verify_bioc.$JOB_ID.out -e logs/verify_bioc.$JOB_ID.err
set -ux
echo "host=$(hostname)"
APPT=$(ls -d /usr/local/package/apptainer/*/bin/apptainer 2>/dev/null|head -1)
export APPTAINER_TMPDIR="${TMPDIR:-/tmp}/apptmp.$JOB_ID"
export APPTAINER_CACHEDIR="${TMPDIR:-/tmp}/apcache.$JOB_ID"
mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"
cd ~/GSE302113
IMG=mgatk_bioc_0.7.0.sif
"$APPT" pull --force "$IMG" docker://quay.io/biocontainers/mgatk:0.7.0--pyhdfd78af_2 && echo "PULL_OK"
echo "=== mgatk version/tenx ==="
"$APPT" exec -B /home "$IMG" mgatk --version 2>&1
"$APPT" exec -B /home "$IMG" mgatk tenx --help >/dev/null 2>&1 && echo "TENX_OK" || echo "TENX_MISSING"
echo "=== samtools ==="
"$APPT" exec -B /home "$IMG" samtools --version 2>&1 | head -1 || echo "SAMTOOLS_MISSING"
echo "=== bedtools ==="
"$APPT" exec -B /home "$IMG" bedtools --version 2>&1 | head -1 || echo "BEDTOOLS_MISSING"
echo "=== R + SummarizedExperiment ==="
"$APPT" exec -B /home "$IMG" R --vanilla -e "suppressMessages({library(SummarizedExperiment);library(Matrix)}); cat(\"R_OK\n\")" 2>&1 | tail -3 || echo "R_MISSING"
echo "=== which Rscript inside ==="
"$APPT" exec -B /home "$IMG" which Rscript R samtools bedtools 2>&1
echo "DONE_VERIFY"
