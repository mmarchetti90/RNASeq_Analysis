process GenerateBowtie2Index {

  label 'bowtie'

  publishDir "${projectDir}/${params.bowtie_index_subdir}", mode: "copy", pattern: "*.bt2"

  input:
  path genome_fasta

  output:
  path "*.bt2", emit: bowtie_index

  """
  bowtie2-build -f ${genome_fasta} ${params.bowtie_index_prefix}
  """

}