#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem=700G
#SBATCH --time=10:00:00
#SBATCH --output=/krummellab/data1/DSCoLab/AMFAR/analysis/10x/2.adt_integration_adtUpDown_rnaUp_cutoffs_v3/adt_integration_adtUpDown_rnaUp_cutoffs_v3_%j.log

set -e
set -o nounset


singularity exec -B /krummellab/ /krummellab/data1/singularity_images/RStudioSingleCell/v4/RStudioSingleCell.sif Rscript /krummellab/data1/DSCoLab/AMFAR/analysis/src/amfar/src/10x/2.adt_integration/adt_integration_adtUpDown_rnaUp_cutoffs_v3.R

