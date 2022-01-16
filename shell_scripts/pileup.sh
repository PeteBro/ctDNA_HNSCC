#!/bin/bash
#SBATCH --job pileup
#SBATCH --ntasks 1
#SBATCH --time 200:0:0
#SBATCH --qos bbdefault
#SBATCH --mem 32G
#SBATCH -o /rds/projects/2017/mehanhhm-group/liquid_biopsies_2021/outputs/variant_calling/job_%A_%a.out
#SBATCH -e /rds/projects/2017/mehanhhm-group/liquid_biopsies_2021/outputs/variant_calling/job_%A_%a.out
set -e

module purge; module load bluebear
module load SAMtools/1.10-GCC-8.3.0
# load modules

workingdir=$(echo "/rds/projects/2017/mehanhhm-group/liquid_biopsies_2021/utilities")
targetdir=$(echo "/rds/projects/2017/mehanhhm-group/liquid_biopsies_2021/data/variant_calling")

echo "$(date) STARTED"

# Running samtools mpileup in a loop over 19 genes
for j in {1..19} ; do
# assuming the bed file has 19 rows for the 19 genes
a=($(sed -n ${j}p ${workingdir}/coords.bed))
# reading the gene coordinates in row number j
reg="${a[0]}:${a[1]}-${a[2]}"
# the region in the form chr:x-y
echo "region $reg"

samtools mpileup --ff 1796 -B -Q0 -q30 -d 10000 -A -f ${workingdir}/hg38.fa -b ${workingdir}/bamlist.txt -r $reg  >> ${targetdir}/allsamples.pileup
	 
done

echo "$(date) FINISHED"