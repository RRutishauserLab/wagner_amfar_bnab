#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=400G
#SBATCH --time=1-00:00:00
#SBATCH --output=/krummellab/data1/DSCoLab/AMFAR/analysis/10x/3.rna_integration_wnnUMAPClust_adtUpDown_rnaUp_cutoffs_v3/rna_integration_wnnUMAPClust_adtUpDown_rnaUp_cutoffs_v3_%j.log


set -e
set -o nounset


singularity exec -B /krummellab/ /krummellab/data1/singularity_images/RStudioSingleCell/v5/RStudioSingleCell.sif Rscript /krummellab/data1/DSCoLab/AMFAR/analysis/src/amfar/src/10x/3.rna_integration_wnnUMAPClust/rna_integration_wnnUMAPClust_adtUpDown_rnaUp_cutoffs_v3.R

