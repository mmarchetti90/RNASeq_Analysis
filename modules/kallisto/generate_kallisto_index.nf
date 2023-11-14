process GenerateKallistoIndex {

  label 'kallisto'

  input:
  path transcripts_fasta

  output:
  path "kallisto_index", emit: kallisto_index

  """
  kallisto index -i kallisto_index ${transcripts_fasta}
  """

}