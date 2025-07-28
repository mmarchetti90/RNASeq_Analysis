#!/usr/bin/env python3

### ---------------------------------------- ###

def load_gtf(path, desired_biotypes=[], desired_chromosomes=[]):
    
    # Load GTF
    
    gtf_data = pd.read_csv(path, sep='\t', header=None, comment='#', dtype=str)
    gtf_data.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

    ### Only keep genes and exons

    gtf_data = gtf_data.loc[gtf_data.feature.isin(['gene', 'exon']), ['seqname', 'start', 'end', 'strand', 'attribute']]
    
    # Get biotype and gene id

    biotypes, gene_ids, gene_symbols, transcript_ids, exon_ids = [], [], [], [], []
    for _,row in gtf_data.iterrows():
        
        info = row.values[-1]
        
        biotype = re.findall('gene_biotype "\w+";', info)[0]
        biotype = biotype.replace('gene_biotype ', '').replace(';', '').replace('"', '')
        
        biotypes.append(biotype)
        
        gene = re.findall('gene_id "\w+";', info)[0]
        gene = gene.replace('gene_id ', '').replace(';', '').replace('"', '')
        
        gene_ids.append(gene)
        
        if 'gene_name' in info:
            
            gene = info[info.index('gene_name "') + len('gene_name "'):]
            gene = gene[:gene.index('"')]
        
        else:
            
            gene = ''
        
        gene_symbols.append(gene)
        
        if 'transcript_id' in info:
            
            transcript = info[info.index('transcript_id "') + len('transcript_id "'):]
            transcript = transcript[:transcript.index('"')]
        
        else:
            
            transcript = ''
        
        transcript_ids.append(transcript)
        
        if 'exon_id' in info:
            
            exon = info[info.index('exon_id "') + len('exon_id "'):]
            exon = exon[:exon.index('"')]
        
        else:
            
            exon = ''
        
        exon_ids.append(exon)

    gtf_data['biotype'] = biotypes
    gtf_data['gene_id'] = gene_ids
    gtf_data['gene_symbol'] = gene_symbols
    gtf_data['transcript_id'] = transcript_ids
    gtf_data['exon_id'] = exon_ids
    
    # Filter based on biotype
    
    if len(desired_biotypes):

        gtf_data = gtf_data.loc[gtf_data.biotype.isin(desired_biotypes),]

    # Filter for desired chromosomes

    if len(desired_chromosomes):

        gtf_data = gtf_data.loc[gtf_data.seqname.isin(desired_chromosomes),]
    
    # Remove genes without gene_symbol
    
    #gtf_data = gtf_data.loc[gtf_data['gene_symbol'] != '',]
    
    # Fix dtypes
    
    gtf_data[['start', 'end']] = gtf_data[['start', 'end']].astype(int)
    
    return gtf_data

### ---------------------------------------- ###

### ------------------MAIN------------------ ###

try:
    
    import numpy as np
    import pandas as pd
    import re
    
    from sys import argv
    
except:
    
    print("One or more dependencies are not installed.\nAlso, make sure your terminal has been activated.")
    exit()

### Load counts

counts_path = argv[argv.index('--counts') + 1]

counts = pd.read_csv(counts_path, sep='\t')

### Load GTF

gtf_path = argv[argv.index('--gtf') + 1]

gtf_data = load_gtf(gtf_path)

gtf_data = gtf_data.loc[(~ gtf_data['exon_id'].isna()) &
                        (gtf_data['exon_id'] != ''),]

gtf_data = gtf_data.drop_duplicates(['seqname', 'start', 'end', 'gene_id', 'exon_id'])

gtf_data.columns = ['Seqname', 'Start', 'End', 'Strand', 'attribute', 'Biotype', 'GeneID', 'GeneSymbol', 'TranscriptID', 'ExonID']

if np.unique(gtf_data['TranscriptID']).shape[0] == 1:
    
    gtf_data = gtf_data.drop(columns=['TranscriptID', 'attribute'])

else:
    
    gtf_data = gtf_data.drop(columns=['attribute'])

### Merge data

annotated_counts = pd.merge(gtf_data, counts, on='ExonID', how='inner')

annotated_counts.to_csv('annotated_counts.tsv', sep='\t', index=False, header=True)
