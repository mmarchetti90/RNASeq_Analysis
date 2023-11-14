process ParseReadsList {

  label 'python'

  publishDir "${params.resources_dir}", mode: "copy", pattern: "*{txt,tsv}"
  
  input:
  path scripts_dir
  path metadata_file

  output:
  path "ReadsList.txt", emit: reads_list
  path "Comparison_*.tsv", optional: true, emit: comparison_files

  """
  python ${scripts_dir}/parse_reads_list/parse_reads_list.py --analysis_type ${workflow.profile} --metadata_file ${metadata_file}
  """

}