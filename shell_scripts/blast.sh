module load engaging/ncbi-blast/2.6.0+

for yfr in 'YFR103' 'YFR107' 'YFR10' 'YFR13' 'YFR2' 'YFR1'
do
    input_file='target_promoter_seqences'/${yfr}'_target_genes_promoter_seqences.fasta'
    output_dir='blast_alignments_output'/${yfr}
    output_file=${output_dir}/${yfr}'_blastn_output.txt'

    echo ${input_file}
    echo ${output_dir}

    blastn -query ${input_file} -subject ${input_file} -out ${output_file}

done