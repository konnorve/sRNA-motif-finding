for yfr in 'YFR103' 'YFR107' 'YFR10' 'YFR13' 'YFR2' 'YFR1'
do
    input_file='target_promoter_BEDs'/${yfr}'_target_genes_promoter_seqences.bed'
    output_dir='gimme_motifs_output'/${yfr}
    genome_ref='genomic_resources/NATL2A_genome_references/onlyNATL2A.fna'

    echo ${input_file}
    echo ${output_dir}

    gimme motifs ${input_file} ${output_dir} --denovo -g ${genome_ref}
done