# RNASeq analysis
## General dockerized Nextflow RNASeq pipeline for bulk and 10X single-cell experiments

/// --------------------------------------- //

### RUN COMMAND:

nextflow run [OPTIONS] --metadata_path "/path/to/metadata/file" rnaseq_analysis.nf

### OPTIONS: (See config file for more)

-profile

		If not specified or equal to "standard", the sample is processed as a bulk RNASeq
		experiment with both a gene-level and transcript-level analyses.

		If equal to "gene_lvl", the sample is processed as a bulk RNASeq experiment with a
		gene-level analysis only.

		If equal to "transcript_lvl", the sample is processed as a bulk RNASeq experiment with a
		transcript-level analysis only.

		If equal to "cellranger", the sample is processed as a single-cell 10X RNASeq experiment
		and UMAPs are produced (models are saved as Pickle files).

		If equal to "small_rna", the sample is processed similalry to "gene_lvl", but with added
		STAR options for better alignment of small-RNA samples.

		If equal to "mirna", the sample is procesed as a miRNA-Seq experiment with Bowtie2
		alignment and annotation of aligned reads to the appropriate mirBase annotation database.

--alignment_only

		If "true", only the alignment step is performed. If "false", the downstream preliminary
		analyses are done as well. Useful if dependencies are not on the server.
		(Default, true)

--trim_reads

		If "true", reads are trimmed with Trim-galore prior to alignment with STAR in standard,
		gene_lvl, and transcript_lvl profiles.
		(Default, true)

--read_length

		Used for STAR index generation.
		(Default, 50)

--saindexnbases

		Used for STAR index generation.
		(Default, 14)

--strand

		For STAR --quantMode counts parsing. Possible values are unstranded, stranded,
		reversestrand.
		(Default, unstranded)
	
--min_counts

		For DESeq2 pre-filtering.
		(Default, 10)

--pval_threshold

		For differentially expressed genes filtering prior to clusterProfiler.
		(Default, 0.05)

--log2fc_threshold

		For differentially expressed genes filtering prior to clusterProfiler.
		(Default, 0.5)

--chemistry

		For Cellranger.
		(Default, auto)

--trimgalore_params

		Additional TrimGalore parameters as a single line of text.
		e.g. "--adapter AAAAAAAAAA --length 15 --clip_R1 3".
		(Default, "")

/// --------------------------------------- ///

### DEPENDENCIES:

Nextflow 20.10+

FastQC

STAR

Kallisto

R 4.1.2+ &

	clusterProfiler
	DESeq2
	enrichplot
	msigdbr
	org.[ORGANISM_OF_INTEREST].eg.db
	pheatmap
	RColorBrewer

Python 3.9.7+ &

	h5py
	igraph
	io
	json
	matplotlib
	numpy
	openpyxl
	os
	pandas
	pickle
	random
	re
	requests
	scipy
	sklearn
	statsmodels
	umap
	seaborn
	shutil
	sys

CellRanger 7+

Bowtie2

Subread featureCounts

Trim-galore &

	Cutadapt

/// --------------------------------------- ///
