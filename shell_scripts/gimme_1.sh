for yfr in 'YFR103' 'YFR107' 'YFR10' 'YFR13' 'YFR2' 'YFR1'
do
    input_file='target_promoter_seqences'/${yfr}'_target_genes_promoter_seqences.fasta'
    output_dir='gimme_motifs_output'/${yfr}

    echo ${input_file}
    echo ${output_dir}

    gimme motifs ${input_file} ${output_dir} --denovo
done