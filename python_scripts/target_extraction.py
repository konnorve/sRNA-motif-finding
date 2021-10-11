from pathlib import Path
from numpy import False_

import pandas as pd
# import numpy as np
import gffpandas.gffpandas as gffpd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

sRNA_targeting_diel_trained_dir = Path('/nobackup1/kve/2021_sRNA_targeting_diel_trained')
genomic_resources_dir = sRNA_targeting_diel_trained_dir / 'genomic_resources'
gff_file = genomic_resources_dir / 'NATL2A_genome_references' / 'onlyNATL2A.gff'
ref_fasta = genomic_resources_dir / 'NATL2A_genome_references' / 'onlyNATL2A.fna'
target_df_path = genomic_resources_dir / 'target_tables' / 'yfrGeneTargets_cleanedUp.csv'

fastas_dir = sRNA_targeting_diel_trained_dir / 'target_promoter_seqences'
beds_dir = sRNA_targeting_diel_trained_dir / 'target_promoter_BEDs'

annotation = gffpd.read_gff3(gff_file)
attributes_df = annotation.attributes_to_columns()

target_df = pd.read_csv(target_df_path)

ref_genome_record = SeqIO.read(ref_fasta, 'fasta')
ref_genome_seq = ref_genome_record.seq

# must match target_df to gff on "new_tag" -- only complete column
# must filter the attrubutes df for cds -- the equivalent column to "new_tag" is "locus_tag"

attributes_df_cds = attributes_df[attributes_df['type']=='CDS']

attributes_df_cds = attributes_df_cds[['seq_id', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase',
                                       'ID', 'Name', 'inference', 'locus_tag', 'product', 'protein_id', 'transl_table']]
attributes_df_cds['new_tag'] = attributes_df_cds['locus_tag']

target_df = target_df.merge(attributes_df_cds, on='new_tag')

target_df['yfr'] = target_df['yfr'].str.replace(' gene targets', '')
target_df['gene_len'] = target_df['end'] - target_df['start']

print(target_df['yfr'].unique())

target_df = target_df.set_index(['yfr', 'new_tag'])

# TODO: practice by pulling out a few protein coding genes and matching up their sequences to the expected values
# TODO: strandedness may be an issue

print(ref_genome_record)
print(type(ref_genome_seq))
print(ref_genome_seq[0])
print(ref_genome_seq[:5])

for yfr in target_df.index.get_level_values('yfr').unique():
    print(yfr)
    yfr_df = target_df.loc[yfr]

    yfr_seq_records = []
    yfr_bed_records = []
    
    for gene in [yfr_df.iloc[x] for x in range(len(yfr_df))]:
        label = gene.name
        start = gene.start
        end = gene.end
        strand = gene.strand

        cds_sequence = ref_genome_seq[start-1:end]
        if strand == '-':
            cds_sequence = cds_sequence.reverse_complement()
        protein_sequence = cds_sequence.translate(table=gene['transl_table'])

        promoter_start = promoter_stop = 0

        if strand == '+':
            promoter_start = start-150
            promoter_stop = start+50
            promoter_region_sequence = ref_genome_seq[promoter_start:promoter_stop]
        elif strand == '-':
            promoter_start = end-49
            promoter_stop = end+151
            promoter_region_sequence = ref_genome_seq[promoter_start:promoter_stop]
            promoter_region_sequence = promoter_region_sequence.reverse_complement()

        print(label, '\t', start, '\t', end, '\t', strand, '\t', cds_sequence[0:5] + '...' + cds_sequence[-5:])

        gene_seq_record = SeqRecord(promoter_region_sequence, 
                                    id=gene.name, name=gene['product'], 
                                    description='promoter region of {} making {} on strand {} that is target of yfr {}'.format(
                                        gene.name, gene['product'], strand, yfr))
        
        # print(gene_seq_record)
        yfr_bed_records.append(['', promoter_start, promoter_stop, label])
        yfr_seq_records.append(gene_seq_record)

    SeqIO.write(yfr_seq_records, fastas_dir / "{}_target_genes_promoter_seqences.fasta".format(yfr), "fasta")
    pd.DataFrame(yfr_bed_records).to_csv(beds_dir / "{}_target_genes_promoter_seqences.bed".format(yfr), sep='\t', index=False, header=False)
