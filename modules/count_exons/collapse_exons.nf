process CollapseExons {
  
  label 'python'

  publishDir "${projectDir}/${params.resources_dir}/${params.collapsed_exons_subdir}", mode: "copy", pattern: "collapsed_exons.gtf"

  input:
  path scripts_dir
  path genome_annotation

  output:
  path "collapsed_exons.gtf", optional: false, emit: collapsed_exons

  """
  ${scripts_dir}/collapse_exons/collapse_exons.py \
  --gtf ${genome_annotation} > collapsed_exons.gtf
  """

}