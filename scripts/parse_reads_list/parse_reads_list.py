#!/usr/bin/env python3

"""
This script reads the metadata file and outputs the following files:
    - List of fastq file names to be passed to alignment (automatically identifies read pairs based on common Sample_IDs);
    - Comparison files for downstream analyses;
"""

### ---------------------------------------- ###

def parse_args():

    # Read analysis type
    if "--analysis_type" not in argv:
    
        sys_exit('ERROR: you need to specify --analysis_type.')

    elif len(argv) <= argv.index("--analysis_type") + 1:
    
        sys_exit('ERROR: you need to specify --analysis_type.')

    else:
    
        analysis_type = argv[argv.index("--analysis_type") + 1]

    # Check that analysis type is valid
    if analysis_type not in ['standard', 'gene_lvl', 'transcript_lvl', 'small_rna', 'mirna', 'cellranger']:
    
        sys_exit('ERROR: invalid --analysis_type.')
    
    # Import project metadata
    if "--metadata_file" not in argv:
        
        sys_exit('ERROR: you need to specify --metadata_file.')
    
    elif len(argv) <= argv.index("--metadata_file") + 1:
        
        sys_exit('ERROR: you need to specify --metadata_file.')
    
    else:
        
        metadata_file_path = argv[argv.index("--metadata_file") + 1]
    
    # Check that metadata file exists
    if not exists(metadata_file_path):
        
        sys_exit(f'ERROR: metadata file {metadata_file_path} not found.')
        
    else:
    
        with open(metadata_file_path) as metadata_file:
    
            metadata = json.load(metadata_file)
    
    # Check that metadata file has correct fields
    if 'sample_data' not in metadata.keys():
        
        sys_exit('ERROR: check that metadata file has the correct fields.')
    
    elif 'sample_id' not in metadata['sample_data'].keys() or 'data_file' not in metadata['sample_data'].keys():
        
        sys_exit(f'ERROR: metadata file {metadata_file_path} not found.')
        
    else:
        
        pass
    
    return analysis_type, metadata

### ---------------------------------------- ###

def listFastqFilesBulk(info):
    
    reads_list = ['SampleID\tFile1\tFile2']
    
    for sample_id in set(info['sample_id']):
        
        read1, *read2 = info.loc[info['sample_id'] == sample_id, 'data_file']
        
        # Manage read2 for paired-end and single-read experiments
        if len(read2):
            
            read2 = read2[0]
            
        else:
            
            read2 = 'mock.fastq'
            with open('mock.fastq', 'w') as mock_file:

                mock_file.write('Empty')
        
        reads_list.append(sample_id + '\t' + read1 + '\t' + read2)
    
    reads_list = '\n'.join(reads_list)

    with open('ReadsList.txt', 'w') as output:

        output.write(reads_list)

### ---------------------------------------- ###

def listFastqFilesSingleCell(info):
    
    reads_list = ['SampleID\tFile1\tFile2']
    
    for sample_id in set(info['sample_id']):
        
        files = list(info.loc[info['sample_id'] == sample_id, 'data_file'])

        if len(files) == 1: # User-provided path is the folder where the fastq files are located
            
            path = files[0]
        
            representative_file = path.split('/')[-1]

        else: # User listed all files (e.g. read1, read2, umi_index)
            
            path = '/'.join(files[0].split('/')[:-1])
            
            representative_file = [f.split('/')[-1] for f in files if '.fastq' in f or 'fq' in f][0]
        
        try:

            sample_prefix = representative_file[:min([search(pattern, representative_file).start() for pattern in ['_S[1-9]*_', '_L[0-9]*_']])]

        except:
            
            sample_prefix = representative_file[:representative_file.index('_')]
        
        reads_list.append(sample_id + '\t' + sample_prefix + '\t' + path)
    
    reads_list = '\n'.join(reads_list)
    
    with open('ReadsList.txt', 'w') as output:

        output.write(reads_list)

### ---------------------------------------- ###

def extrapolateGroups(sample_ids):
    
    """
    This function will try to guess sample groups assuming that samples from the same group have the same naming convention
    """
    
    # Splitting Sample_IDs into words (separated by either spaces or underscores)
    groups = [i.replace(' ', '_').split('_') for i in sample_ids]
    
    # Removing isolated numbers in range 0-1000 (which most likely indicate condition replicate number), then re-joining words
    string_numbers = [str(n) for n in range(1001)]
    
    for i,g in enumerate(groups):
        
        parsed_id = '_'.join([word for word in g if word not in string_numbers])
        groups[i] = parsed_id
    
    # Adding 'Group' column
    updated_info = pd.DataFrame({'Sample_ID' : sample_ids, 'Group' : groups})
    
    # Removing groups made of only one element
    group_filter = [groups.count(g) > 1 for g in groups]
    updated_info = updated_info.loc[group_filter,]
    
    return updated_info

### ---------------------------------------- ###

def createComparisonFiles(info):
    
    keywords = ['control', 'con', 'ctrl', 'cneg', 'c-neg', 'wildtype', 'wild-type', 'wt'] # Possible control sample names
    
    comparisons = []
    for g1 in info['Group'].unique():

        for g2 in info['Group'].unique():
            
            # mMking sure the comparison is new and that g1 != g2
            if g1 + '_vs_' + g2 in comparisons or g1 == g2:

                continue
            
            comparisons.append(g1 + '_vs_' + g2)
            comparisons.append(g2 + '_vs_' + g1)
            
            # Determining reference by looking for keywords in the g1 group name
            if sum([key in g1.lower() for key in keywords]) != 0:
                
                reference = g1
                comparison_name = g2 + '_vs_' + g1
                
            else:
                
                reference = g2
                comparison_name = g1 + '_vs_' + g2
            
            # Creating new comparison file
            new_comparison = ['ComparisonName\t' + comparison_name,
                              'ComparisonDesing\t~ condition',
                              'Reference\t' + reference,
                              'sample\tcondition']
            new_comparison += set(['\t'.join(row) for _,row in info.loc[(info.Group == g1) | (info.Group == g2)].iterrows()])
            new_comparison = '\n'.join(new_comparison)
            
            # Saving to file
            with open('Comparison_' + comparison_name + '.tsv', 'w') as output:

                output.write(new_comparison)

### ------------------MAIN------------------ ###

import json
import pandas as pd

from os.path import exists
from re import search
from sys import argv
from sys import exit as sys_exit

# Parse arguments
analysis_type, metadata = parse_args()

# Read sample table
sample_data = pd.DataFrame(metadata['sample_data'])

# Making a list of fastq files to be aligned/counted
if analysis_type in ['standard', 'gene_lvl', 'transcript_lvl', 'small_rna', 'mirna']:

    listFastqFilesBulk(sample_data)

else:
    
    listFastqFilesSingleCell(sample_data)

# Creating DESeq2 comparisons files for gene-level and small-rna analysis
if analysis_type in ['standard', 'gene_lvl', 'transcript_lvl', 'small_rna', 'mirna']:

    groups = extrapolateGroups(sample_data.sample_id.unique())

    createComparisonFiles(groups)
