#!/bin/bash
#SBATCH --job vcf
#SBATCH --ntasks 1
#SBATCH --time 48:0:0
#SBATCH --qos bbdefault
#SBATCH --mem 32G
#SBATCH -o /rds/projects/2017/mehanhhm-group/liquid_biopsies_2021/outputs/variant_calling/job_%A_%a.out
#SBATCH -e /rds/projects/2017/mehanhhm-group/liquid_biopsies_2021/outputs/variant_calling/job_%A_%a.out
set -e


module purge; module load bluebear
module load VarScan/2.4.4-Java-1.8
# load modules

workingdir=$(echo "/rds/projects/2017/mehanhhm-group/liquid_biopsies_2021/data/variant_calling")

java -jar $EBROOTVARSCAN/VarScan.v2.4.4.jar mpileup2cns ${workingdir}/allsamples.pileup --vcf-sample-list samples.txt --output-vcf 1 --strand-filter 0 --min-avg-qual 0 --variants 1 --min-coverage 0 --min-reads2 0 --min-var-freq 0.01 --p-value 1 1> ${workingdir}/allsamples 2> ${workingdir}/allsamples.log
# search for variants