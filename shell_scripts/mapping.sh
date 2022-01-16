#!/bin/bash
#SBATCH --job mapping
#SBATCH --ntasks 1
#SBATCH --time 24:0:0
#SBATCH --qos bbdefault
#SBATCH --mem 32G
#SBATCH --array 1-50
#SBATCH -o /rds/projects/2017/mehanhhm-group/liquid_biopsies_2021/outputs/mapping/logs/job_%A_%a.out
#SBATCH -e /rds/projects/2017/mehanhhm-group/liquid_biopsies_2021/outputs/mapping/logs/job_%A_%a.out
BB_WORKDIR=$(mktemp -d /scratch/${USER}_${SLURM_JOBID}.XXXXXX)
export TMPDIR=${BB_WORKDIR}
set -e

module purge; module load bluebear
module load BWA/0.7.17-GCC-8.3.0
module load samblaster/0.1.26-GCCcore-8.3.0
module load SAMtools/1.10-GCC-8.3.0
# load modules

workingdir=$(echo "/rds/projects/2017/mehanhhm-group/liquid_biopsies_2021/data/trimmed_fastq")
targetdir=$(echo "/rds/projects/2017/mehanhhm-group/liquid_biopsies_2021/data/mapped_fastq")
reference_genome=$(echo "/rds/projects/2017/mehanhhm-group/liquid_biopsies_2021/utilities/hg38.fa")
sampledir=$(echo "/rds/projects/2017/mehanhhm-group/liquid_biopsies_2021/utilities")
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${sampledir}/samples.txt)
bwa mem ${reference_genome} ${workingdir}/*${sample}_R1* ${workingdir}/*${sample}_R2* | samblaster | samtools view -bh | samtools sort > ${targetdir}/${sample}.bam
# align reads to reference genome with bwa
# mark duplicates with samblaster

test -d ${BB_WORKDIR} && /bin/rm -rf ${BB_WORKDIR}