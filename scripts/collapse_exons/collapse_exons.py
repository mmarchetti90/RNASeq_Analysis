#!/usr/bin/env python3

"""
Collapse exons for DEXseq
"""

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

def collapse_exons(dt):
    
    bad_junction = 'end_start'
    
    for seqname in np.unique(dt['seqname']):
        
        for strand in ['+', '-']:
            
            # Subset dt
            
            dt_sub = dt.loc[(dt['seqname'] == seqname) & (dt['strand'] == strand),]
            
            # Create a dataframe where start and end positions of exons are sorted together
            
            pos_dt = pd.DataFrame({c : [] for c in ['gene_id', 'gene_symbol', 'biotype']})
            
            for col in ['start', 'end']:
                
                pos_sub = dt_sub[['gene_id', 'gene_symbol', 'biotype']].copy()
                pos_sub.loc[:, 'pos'] = dt_sub[col].values
                pos_sub.loc[:, 'pos_type'] = [col] * pos_sub.shape[0]
                
                pos_dt = pd.concat([pos_dt, pos_sub], axis=0, ignore_index=True)
            
            # Sort by coordinate
            
            pos_dt = pos_dt.sort_values('pos').reset_index(drop=True)
            
            #pos_sub.loc[:, 'pos'] = pos_sub['pos'].astype(int)
            
            # Add entries to GTF
            
            """
            # The following has issues when one gene is "inside" another.
            # e.g. ENSG00000000971 and ENSG00000289697
            
            exon_count = 1
            
            previous_gene_id = ''
            
            for (_,row1),(_,row2) in zip(pos_dt.iterrows(), pos_dt.iloc[1:,].iterrows()):
                
                bin_type = row1['pos_type'] + '_' + row2['pos_type']
                
                if bin_type != bad_junction and row1['gene_id'] == row2['gene_id']:
                    
                    start, end = str(int(row1['pos'])), str(int(row2['pos']))
                    
                    gene_id, gene_symbol, gene_biotype = row1[['gene_id', 'gene_symbol', 'biotype']].values
                    
                    if gene_id == previous_gene_id:
                        
                        exon_count += 1
                    
                    else:
                        
                        exon_count = 1
                        
                        previous_gene_id = gene_id
            
                    if gene_symbol != '':
                        
                        attribute = f'gene_id "{gene_id}"; gene_name "{gene_symbol}"; gene_biotype "{gene_biotype}"; exon_number "{exon_count}"; exon_id "{gene_id}_E{exon_count}"; exon_bin_type "{bin_type}"'
                    
                    else:
                        
                        attribute = f'gene_id "{gene_id}"; gene_biotype "{gene_biotype}"; exon_number "{exon_count}"; exon_id "{gene_id}_E{exon_count}"; exon_bin_type "{bin_type}"'
                    
                    print('\t'.join([seqname, 'collapse_exons', 'exon', start, end, '.', strand, '.', attribute]))
            
            # The following code, processing 1 gene at a time, solves this
            """
            
            for gene_id in np.unique(pos_dt['gene_id']):
                
                pos_dt_sub = pos_dt.loc[pos_dt['gene_id'] == gene_id]
                
                gene_id, gene_symbol, gene_biotype = pos_dt_sub.iloc[0][['gene_id', 'gene_symbol', 'biotype']].values
                
                exon_count = 0
                
                for (_,row1),(_,row2) in zip(pos_dt_sub.iterrows(), pos_dt_sub.iloc[1:,].iterrows()):
                    
                    bin_type = row1['pos_type'] + '_' + row2['pos_type']
                    
                    if bin_type != bad_junction and row1['gene_id'] == row2['gene_id']:
                        
                        exon_count += 1
                        
                        start, end = str(int(row1['pos'])), str(int(row2['pos']))
                        
                        if gene_symbol != '':
                            
                            attribute = f'gene_id "{gene_id}"; gene_name "{gene_symbol}"; gene_biotype "{gene_biotype}"; exon_number "{exon_count}"; exon_id "{gene_id}_E{exon_count}"; exon_bin_type "{bin_type}"'
                        
                        else:
                            
                            attribute = f'gene_id "{gene_id}"; gene_biotype "{gene_biotype}"; exon_number "{exon_count}"; exon_id "{gene_id}_E{exon_count}"; exon_bin_type "{bin_type}"'
                        
                        print('\t'.join([seqname, 'collapse_exons', 'exon', start, end, '.', strand, '.', attribute]))

### ------------------MAIN------------------ ###

try:
    
    import numpy as np
    import pandas as pd
    import re
    
    from sys import argv
    
except:
    
    print("One or more dependencies are not installed.\nAlso, make sure your terminal has been activated.")
    exit()

### Load GTF

gtf_path = argv[argv.index('--gtf') + 1]

gtf_data = load_gtf(gtf_path)

gtf_data = gtf_data.loc[(~ gtf_data['exon_id'].isna()) &
                        (gtf_data['exon_id'] != ''),]

gtf_data = gtf_data.drop_duplicates(['seqname', 'start', 'end', 'gene_id'])

### Collapse exons

collapsed_gtf_data = collapse_exons(gtf_data)
