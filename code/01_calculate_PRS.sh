#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# POLYGENIC RISK SCORE CALCULATION
#
# Manuscript:
# Polygenic Risk Scores for Hypertensive Disorders of
# Pregnancy and Later-Life Cardiovascular Disease in Women
#
# Software:
#   pgsc_calc v2.2.0
#   Nextflow v25.04.6
#   Singularity v4.1.0
#
# Target genome build:
#   GRCh38
#
# Scores:
#   PGS003586
#   PGS003587
#   PGS004593
#   Broad HDP PRS
#
# Replace the example paths below with paths appropriate
# for the local computing environment.
# ============================================================


# ============================================================
# PATHS
# ============================================================

SAMPLESHEET_PE_GH="/path/to/samplesheet_PE_GH.csv"

SAMPLESHEET_HDP="/path/to/samplesheet_HDP.csv"

HDP_SCORE_FILE="/path/to/015_HYPTENSPREG_FinnGen_r8_PRScs_auto.pgsc_calc.txt"

ANCESTRY_REFERENCE="/path/to/pgsc_HGDP+1kGP_v1.tar.zst"

HG19_TO_HG38="/path/to/hg19ToHg38.over.chain.gz"

HG38_TO_HG19="/path/to/hg38ToHg19.over.chain.gz"

OUTDIR="/path/to/results"

WORKDIR="/path/to/nextflow_work"


# ============================================================
# 1. PGS CATALOG SCORES
#
# PGS003586
# PGS003587
# PGS004593
# ============================================================

nextflow run pgscatalog/pgsc_calc \
    -r v2.2.0 \
    -profile singularity \
    --input "${SAMPLESHEET_PE_GH}" \
    --pgs_id PGS003586,PGS003587,PGS004593 \
    --target_build GRCh38 \
    --run_ancestry "${ANCESTRY_REFERENCE}" \
    --outdir "${OUTDIR}" \
    -work-dir "${WORKDIR}"


# ============================================================
# 2. BROAD HDP PRS
#
# FinnGen-based HDP score.
#
# LiftOver was performed within pgsc_calc before scoring
# against the GRCh38 target genotype data.
# ============================================================

nextflow run pgscatalog/pgsc_calc \
    -r v2.2.0 \
    -profile singularity \
    --input "${SAMPLESHEET_HDP}" \
    --scorefile "${HDP_SCORE_FILE}" \
    --liftover \
    --target_build GRCh38 \
    --hg19_chain "${HG19_TO_HG38}" \
    --hg38_chain "${HG38_TO_HG19}" \
    --run_ancestry "${ANCESTRY_REFERENCE}" \
    --outdir "${OUTDIR}" \
    -work-dir "${WORKDIR}"