#!/usr/bin/env bash
# 
#SBATCH --job-name=gimme_motifs
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
#SBATCH -o logs/%j_%a_slurm_output.txt
#SBATCH -e logs/%j_%a_slurm_error.txt

conda activate gimme

working_dir="/nobackup1/kve/2021_sRNA_targeting_diel_trained"

for yfr in 'YFR103' 'YFR107' 'YFR10' 'YFR13' 'YFR2' 'YFR1'
do
    input_file=${working_dir}/'target_promoter_seqences'/${yfr}'_target_genes_promoter_seqences.fasta'
    output_dir=${working_dir}/'gimme_motifs_output'/${yfr}
    genome_ref=${working_dir}/'genomic_resources/NATL2A_genome_references/onlyNATL2A.fna'

    echo ${input_file}
    echo ${output_dir}

    gimme motifs ${input_file} ${output_dir} --denovo -g ${genome_ref}
done