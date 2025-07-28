process GenerateKallistoIndex {

  label 'kallisto'

  publishDir "${projectDir}/${params.resources_dir}/${params.kallisto_index_subdir}", mode: "copy", pattern: "kallisto_index"

  input:
  path transcripts_fasta

  output:
  path "kallisto_index", emit: kallisto_index

  """
  kallisto index -t \$SLURM_CPUS_ON_NODE -i kallisto_index ${transcripts_fasta}
  """

}