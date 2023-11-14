process CountExons {
  
  label 'subread'

  publishDir "${projectDir}/${params.reports_dir}/${params.exons_summary}", mode: "copy", pattern: "*.AlignmentSummary.tsv"
  publishDir "${projectDir}/${params.exon_counts_dir}", mode: "copy", pattern: "*_exon_counts.tsv"

  input:
  each path(genome_annotation)
  tuple val(read_id), path(bam), val(read2_name)

  output:
  path "${read_id}_exon_counts.tsv", optional: true, emit: exon_counts
  path "${read_id}.AlignmentSummary.tsv", optional: true, emit: exon_counts_report

  """
  # Strand information
  if [[ "${params.strand}" == "stranded" ]]
  then

      strand_type=1

  elif [[ "${params.strand}" == "reversestrand" ]]
  then

      strand_type=2

  else

      strand_type=0

  fi

  # Running featureCounts
  if [[ "${read2_name}" == "mock.fastq" ]]
  then

      # Counting exons with Subread featureCounts for single-end reads
      featureCounts -f -t exon -g exon_id -s \$strand_type -O -a ${genome_annotation} -F GTF -o temp.tsv ${bam}

  else

      # Counting exons with Subread featureCounts for paired-end reads
      featureCounts -p -f -t exon -g exon_id -s \$strand_type -O -a ${genome_annotation} -F GTF -o temp.tsv ${bam}

  fi

  # Removing first line of count file (it contains the command used to run featureCounts)
  sed "1d" temp.tsv > ${read_id}_exon_counts.tsv

  # Renaming summary from featureCounts
  mv temp.tsv.summary ${read_id}.AlignmentSummary.tsv
  """

}