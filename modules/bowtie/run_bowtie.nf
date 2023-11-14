process RunBowtie2 {

  // Using Bowtie2 to align miRNA (works best according to Ziemann et al., "Evaluation of microRNA alignment techniques", RNA, 2016)
  // Reads are aligned to the whole genome (aligning only to mirBase mature.fa list of miRNA sequences does not seem to work well)
  // Aligned reads are then counted with featureCounts using a mirBase annotation

  label 'bowtie'

  publishDir "${projectDir}/${params.bam_dir}", mode: "copy", pattern: "*.bam"
  publishDir "${projectDir}/${params.reports}/${params.bowtie_summary}", mode: "copy", pattern: "*.AlignmentSummary.tsv"
  publishDir "${projectDir}/${params.gene_counts_dir}", mode: "copy", pattern: "*_miRNA_counts.tsv"

  input:
  path(bowtie_index)
  tuple val(read_id), path(read1), path(read2)
  each path(mirna_annotation)

  output:
  path "${read_id}_Aligned.sortedByCoord.out.bam", optional: true
  path "${read_id}_miRNA_counts.tsv", optional: true, emit: mirna_counts
  path "*.AlignmentSummary.tsv", emit: alignment_summary

  """
  # Bowtie alignment on the whole genome
  if [[  "${read2})" == "mock.fastq" ]]
  then

    bowtie2 --threads \$SLURM_CPUS_ON_NODE --local --very-sensitive-local -x bowtie_index -U ${read1} -S ${read_id}_Aligned.out.sam

  else

    bowtie2 --threads \$SLURM_CPUS_ON_NODE --local --very-sensitive-local -x bowtie_index -1 ${read1} -2 ${read2} -S ${read_id}_Aligned.out.sam

  fi

  # Sorting by coordinates
  samtools sort ${read_id}_Aligned.out.sam -o ${read_id}_Aligned.sortedByCoord.out.sam
  
  # Convert sam to bam
  samtools view --threads \$SLURM_CPUS_ON_NODE -S -b ${read_id}_Aligned.sortedByCoord.out.sam > ${read_id}_Aligned.sortedByCoord.out.bam

  # Counting miRNAs with Subread featureCounts
  featureCounts -t miRNA -g Name -O -s 1 -a ${mirna_annotation} -F GTF -o temp.tsv ${read_id}_Aligned.sortedByCoord.out.sam

  # Removing first line of count file (it contains the command above)
  sed "1d" temp.tsv > ${read_id}_miRNA_counts.tsv

  # Renaming summary from featureCounts
  mv temp.tsv.summary ${read_id}.AlignmentSummary.tsv
  """

}