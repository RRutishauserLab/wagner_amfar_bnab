#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --time=4:00:00
#SBATCH --output=/krummellab/data1/DSCoLab/AMFAR/analysis/10x/1.dsb_normalization_adtUpDown_rnaUp_cutoffs_v2/dsb_normalization_%j.log

set -e
set -o nounset


singularity exec -B /krummellab/ /krummellab/data1/singularity_images/RStudioSingleCell/v4/RStudioSingleCell.sif Rscript /krummellab/data1/DSCoLab/AMFAR/analysis/src/amfar/src/10x/1.dsb_normalization/perform_dsb_normalization_adtUpDown_rnaUp_cutoffs_v2.R $SLURM_ARRAY_TASK_ID

