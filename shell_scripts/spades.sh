module load engaging/SPAdes/3.6.2
module load engaging/python/3.5.1

for yfr in 'YFR103' 'YFR107' 'YFR10' 'YFR13' 'YFR2' 'YFR1'
do
    input_file='target_promoter_seqences'/${yfr}'_target_genes_promoter_seqences.fasta'
    output_dir='spades_alignments'/${yfr}

    echo ${input_file}
    echo ${output_dir}

    spades.py -s ${input_file} -o ${output_dir} --only-assembler --cov 1
done

module unload engaging/SPAdes/3.6.2
module unload engaging/python/3.5.1