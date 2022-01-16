#!/bin/bash
#SBATCH --job annotate
#SBATCH --ntasks 1
#SBATCH --time 24:0:0
#SBATCH --qos bbdefault
#SBATCH --mem 16G
#SBATCH -o /rds/projects/2017/mehanhhm-group/liquid_biopsies_2021/outputs/variant_calling/job_%A_%a.out
#SBATCH -e /rds/projects/2017/mehanhhm-group/liquid_biopsies_2021/outputs/variant_calling/job_%A_%a.out
set -e

module purge; module load bluebear
# load modules

workingdir=$(echo "/rds/projects/2017/mehanhhm-group/liquid_biopsies_2021/utilities")
targetdir=$(echo "/rds/projects/2017/mehanhhm-group/liquid_biopsies_2021/data/variant_calling")

perl $workingdir/annovar/annovar/table_annovar.pl $targetdir/allsamples.vcf $workingdir/annovar/annovar/humandb/ --vcfinput -v --nastring . --buildver hg38 --protocol cosmic94_noncoding,cosmic94_coding -operation f,f --remove
