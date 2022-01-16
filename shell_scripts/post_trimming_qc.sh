#!/bin/bash
#SBATCH --job QC
#SBATCH --ntasks 1
#SBATCH --time 0:10:0
#SBATCH --qos bbshort
#SBATCH --mem 8G
#SBATCH --array 1-100
#SBATCH -o /rds/projects/2017/mehanhhm-group/liquid_biopsies_2021/outputs/post_trimming_qc/logs/job_%A_%a.out
#SBATCH -e /rds/projects/2017/mehanhhm-group/liquid_biopsies_2021/outputs/post_trimming_qc/logs/job_%A_%a.out
set -e

module purge; module load bluebear
module load FastQC/0.11.9-Java-11
# load modules

workingdir=$(echo "/rds/projects/2017/mehanhhm-group/liquid_biopsies_2021/data/trimmed_fastq")
filename=$(ls ${workingdir} | sed -n ${SLURM_ARRAY_TASK_ID}p)
targetdir=$(echo "/rds/projects/2017/mehanhhm-group/liquid_biopsies_2021/outputs/post_trimming_qc/$(echo ${filename} | cut -d"." -f1)")
mkdir ${targetdir}
echo "job ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_TASK_MAX}: ${filename}"
fastqc ${workingdir}/${filename} -o ${targetdir}
# run fastqc quality-control check on all sequence data files