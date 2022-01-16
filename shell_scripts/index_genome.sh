#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --time 6:0:0
#SBATCH --qos bbdefault
#SBATCH --mem 8G
set -e

module purge; module load bluebear
module load BWA/0.7.17-foss-2019a
# load modules

bwa index -a bwtsw /rds/projects/2017/mehanhhm-group/liquid_biopsies_2021/utilities/hg38.fa
# index reference genome