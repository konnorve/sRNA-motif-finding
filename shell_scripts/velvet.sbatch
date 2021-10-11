#!/usr/bin/env bash
# 
#SBATCH --job-name=velvet
#
# Specifies using a centos7 node
#SBATCH -C centos7
#
# wall clock limit:
#SBATCH --time 24:00:00
#
# Partition to use:
#SBATCH --partition sched_mit_chisholm
#
# Number of tasks/cores for job
#SBATCH -n 1
#
#SBATCH --array=1,2,13,10,103,107
#
#SBATCH --mem 16000
#
#SBATCH -o logs/velvet/%j_%a_slurm_output.txt
#SBATCH -e logs/velvet/%j_%a_slurm_error.txt

# 'YFR103' 'YFR107' 'YFR10' 'YFR13' 'YFR2' 'YFR1'
YFR="YFR${SLURM_ARRAY_TASK_ID}"

KMER_LEN=7

working_dir="/nobackup1/kve/2021_sRNA_targeting_diel_trained"

input_file=${working_dir}/'target_promoter_seqences'/${YFR}'_target_genes_promoter_seqences.fasta'
output_dir=${working_dir}/'velvet_output'/${YFR}
genome_ref=${working_dir}/'genomic_resources/NATL2A_genome_references/onlyNATL2A.fna'

echo ${input_file}
echo ${output_dir}

mkdir -p ${output_dir}

velveth ${output_dir} ${KMER_LEN} -short -fasta ${input_file} 

velvetg ${output_dir}