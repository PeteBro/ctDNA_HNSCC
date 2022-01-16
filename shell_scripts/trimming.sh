#!/bin/bash
#SBATCH --job trimming
#SBATCH --ntasks 1
#SBATCH --time 24:0:0
#SBATCH --qos bbdefault
#SBATCH --mem 16G
#SBATCH --array 1-50
#SBATCH -o /rds/projects/2017/mehanhhm-group/liquid_biopsies_2021/outputs/trimming/job_%A_%a.out
#SBATCH -e /rds/projects/2017/mehanhhm-group/liquid_biopsies_2021/outputs/trimming/job_%A_%a.out
BB_WORKDIR=$(mktemp -d /scratch/${USER}_${SLURM_JOBID}.XXXXXX)
export TMPDIR=${BB_WORKDIR}
set -e

module purge; module load bluebear
module load Trimmomatic/0.39-Java-11
# load modules

workingdir=$(echo "/rds/projects/2017/mehanhhm-group/liquid_biopsies_2021/data/raw_fastq")
targetdir=$(echo "/rds/projects/2017/mehanhhm-group/liquid_biopsies_2021/data/trimmed_fastq")
samples=$(echo "/rds/projects/2017/mehanhhm-group/liquid_biopsies_2021/utilities/samples.txt")
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${samples})
echo "job ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_TASK_MAX}: ${sample}"
echo "input files: $(ls ${workingdir}/*${sample}_R1* | rev | cut -d "/" -f 1 | rev), $(ls ${workingdir}/*${sample}_R2* | rev | cut -d "/" -f 1 | rev)"

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -phred33 ${workingdir}/*${sample}_R1* ${workingdir}/*${sample}_R2* ${targetdir}/${sample}_R1P.fastq.gz ${targetdir}/${sample}_R1U.fastq.gz ${targetdir}/${sample}_R2P.fastq.gz ${targetdir}/${sample}_R2U.fastq.gz LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:30
# run paired end trimming on reads corresponding to each sample

test -d ${BB_WORKDIR} && /bin/rm -rf ${BB_WORKDIR}