process RunSTAR {

  // Similar to the same process in gene_lvl workflow, but with some additional parameters that according to STAR creator should help with small RNAs
  
  label 'star'

  publishDir "${projectDir}/${params.bam_dir}", mode: "copy", pattern: "*{Aligned,Unmapped,SJ}*"
  publishDir "${projectDir}/${params.gene_counts_dir}", mode: "copy", pattern: "*_ReadsPerGene.out.tab"
  publishDir "${projectDir}/${params.reports_dir}/${params.star_reports}", mode: "copy", pattern: "*Log*"

  input:
  each path(index)
  tuple val($read_id), val(read_id), path(read1), path(read2)

  output:
  path "${read_id}_Aligned.sortedByCoord.out.bam", optional: true, emit: bam_files
  path "${read_id}_Unmapped.out.mate1", optional: true
  path "${read_id}_ReadsPerGene.out.tab", optional: true, emit: gene_counts
  path "${read_id}_Log.{final.out,out,progress.out}", optional: true, emit: star_reports
  path "${read_id}_SJ.out.tab", optional: true
  
  """
  # STAR alignment
  temp="${read1}"

  if [[ \${temp: -2} == "gz" ]]
  then

    if [[ "${read2}" == "mock.fastq" ]]
    then

      STAR \
      --runThreadN \$SLURM_CPUS_ON_NODE \
      --runMode alignReads \
      --twopassMode Basic \
      --genomeDir ${index} \
      --readFilesIn ${read1} \
      --outFileNamePrefix ./${read_id}_ \
      --outSAMtype BAM Unsorted \
      --outReadsUnmapped Fastx \
      --quantMode GeneCounts \
      --readFilesCommand zcat \
      --alignEndsType EndToEnd \
      --alignIntronMax 1 \
      --outFilterMatchNmin 16 \
      --alignSJDBoverhangMin 10000 \
      --outFilterMismatchNoverLmax 0.05 \
      --outFilterScoreMinOverLread 0 \
      --outFilterMatchNminOverLread 0

    else

      STAR \
      --runThreadN \$SLURM_CPUS_ON_NODE \
      --runMode alignReads \
      --twopassMode Basic \
      --genomeDir ${index} \
      --readFilesIn ${read1} ${read2} \
      --outFileNamePrefix ./${read_id}_ \
      --outSAMtype BAM Unsorted \
      --outReadsUnmapped Fastx \
      --quantMode GeneCounts \
      --readFilesCommand zcat \
      --alignEndsType EndToEnd \
      --alignIntronMax 1 \
      --outFilterMatchNmin 16 \
      --alignSJDBoverhangMin 10000 \
      --outFilterMismatchNoverLmax 0.05 \
      --outFilterScoreMinOverLread 0 \
      --outFilterMatchNminOverLread 0

    fi

  else

    if [[ "${read2}" == "mock.fastq" ]]
    then

      STAR \
      --runThreadN \$SLURM_CPUS_ON_NODE \
      --runMode alignReads \
      --twopassMode Basic \
      --genomeDir ${index} \
      --readFilesIn ${read1} \
      --outFileNamePrefix ./${read_id}_ \
      --outSAMtype BAM Unsorted \
      --outReadsUnmapped Fastx \
      --quantMode GeneCounts \
      --alignEndsType EndToEnd \
      --alignIntronMax 1 \
      --outFilterMatchNmin 16 \
      --alignSJDBoverhangMin 10000 \
      --outFilterMismatchNoverLmax 0.05 \
      --outFilterScoreMinOverLread 0 \
      --outFilterMatchNminOverLread 0

    else

      STAR \
      --runThreadN \$SLURM_CPUS_ON_NODE \
      --runMode alignReads \
      --twopassMode Basic \
      --genomeDir ${index} \
      --readFilesIn ${read1} ${read2} \
      --outFileNamePrefix ./${read_id}_ \
      --outSAMtype BAM Unsorted \
      --outReadsUnmapped Fastx \
      --quantMode GeneCounts \
      --alignEndsType EndToEnd \
      --alignIntronMax 1 \
      --outFilterMatchNmin 16 \
      --alignSJDBoverhangMin 10000 \
      --outFilterMismatchNoverLmax 0.05 \
      --outFilterScoreMinOverLread 0 \
      --outFilterMatchNminOverLread 0

    fi

  fi

  # Sorting by coordinates
  # N.B. I use samtools instead of STAR --outSAMtype BAM SortedByCoordinate because STAR was sometimes runnin into memory issues on the cluster
  samtools sort ${read_id}_Aligned.out.bam -o ${read_id}_Aligned.sortedByCoord.out.bam -@ \$SLURM_CPUS_ON_NODE
  """

}