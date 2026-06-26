#!/bin/bash
## Submit all ALS-iMN-ML analysis jobs to LSF (Compute1)
## Usage: bash submit/submit_all.sh

REPO_DIR="/rdcw/fs1/jmilbrandt/Active/Neuronal_Resilience_Program/ALS-iMN-ML"
BRBSEQ_DIR="/rdcw/fs1/jmilbrandt/Active/Neuronal_Resilience_Program/NRP-BRBseq"
STORE_DIR="/storage1/fs1/jmilbrandt/Active/Neuronal_Resilience_Program"

export LSF_DOCKER_VOLUMES="\
${REPO_DIR}:${REPO_DIR} \
${BRBSEQ_DIR}:${BRBSEQ_DIR} \
${STORE_DIR}:${STORE_DIR}"

DOCKER="docker(bioconductor/bioconductor_docker:RELEASE_3_20)"
LOG_DIR="${REPO_DIR}/logs"
mkdir -p "${LOG_DIR}"

submit_job() {
  local JOB_NAME=$1
  local MEM=$2
  local SCRIPT=$3
  bsub \
    -J  "${JOB_NAME}" \
    -G  compute-jmilbrandt \
    -q  general \
    -n  4 \
    -M  "${MEM}" \
    -R  "rusage[mem=${MEM}] span[hosts=1]" \
    -W  1:00 \
    -o  "${LOG_DIR}/${JOB_NAME}_%J.log" \
    -a  "${DOCKER}" \
    -env "REPO_DIR=${REPO_DIR}" \
    Rscript "${REPO_DIR}/scripts/${SCRIPT}"
}

# Job 1 must complete before Jobs 2-6 (it generates feature RDS files)
JOB1=$(submit_job "als_prepare"    16GB "01_prepare_features.R")
echo "Submitted: ${JOB1}"
JOB1_ID=$(echo "${JOB1}" | grep -oP '(?<=Job <)\d+')

# Jobs 2-4 depend on features (parallel)
bsub -J "als_enet_binary"  -G compute-jmilbrandt -q general -n 4 -M 32GB \
  -R "rusage[mem=32GB] span[hosts=1]" -W 1:30 \
  -w "done(${JOB1_ID})" \
  -o "${LOG_DIR}/als_enet_binary_%J.log" \
  -a "${DOCKER}" -env "REPO_DIR=${REPO_DIR}" \
  Rscript "${REPO_DIR}/scripts/02_elastic_net_binary.R"

bsub -J "als_reg_compare"  -G compute-jmilbrandt -q general -n 4 -M 32GB \
  -R "rusage[mem=32GB] span[hosts=1]" -W 1:30 \
  -w "done(${JOB1_ID})" \
  -o "${LOG_DIR}/als_reg_compare_%J.log" \
  -a "${DOCKER}" -env "REPO_DIR=${REPO_DIR}" \
  Rscript "${REPO_DIR}/scripts/03_regularization_compare.R"

bsub -J "als_rf"           -G compute-jmilbrandt -q general -n 4 -M 32GB \
  -R "rusage[mem=32GB] span[hosts=1]" -W 1:30 \
  -w "done(${JOB1_ID})" \
  -o "${LOG_DIR}/als_rf_%J.log" \
  -a "${DOCKER}" -env "REPO_DIR=${REPO_DIR}" \
  Rscript "${REPO_DIR}/scripts/04_random_forest.R"

bsub -J "als_de_enet"      -G compute-jmilbrandt -q general -n 4 -M 32GB \
  -R "rusage[mem=32GB] span[hosts=1]" -W 1:30 \
  -w "done(${JOB1_ID})" \
  -o "${LOG_DIR}/als_de_enet_%J.log" \
  -a "${DOCKER}" -env "REPO_DIR=${REPO_DIR}" \
  Rscript "${REPO_DIR}/scripts/05_de_enet_pipeline.R"

bsub -J "als_temporal"     -G compute-jmilbrandt -q general -n 4 -M 32GB \
  -R "rusage[mem=32GB] span[hosts=1]" -W 1:30 \
  -w "done(${JOB1_ID})" \
  -o "${LOG_DIR}/als_temporal_%J.log" \
  -a "${DOCKER}" -env "REPO_DIR=${REPO_DIR}" \
  Rscript "${REPO_DIR}/scripts/06_temporal_analysis.R"

echo ""
echo "All jobs submitted. Monitor with: bjobs"
echo "Logs -> ${LOG_DIR}/"
