process IntegrateAnalyses {

  label 'r'

  publishDir "${projectDir}/${params.dexseq_sleuth_integration}/${params.common_genes_dir}", mode: 'copy', pattern: "CommonGenes_*.tsv"

  input:
  each path(scripts_dir)
  each path(genome_annotation)
  tuple val(analysis_name), path(dexseq_file), path(sleuth_file)

  output:
  path "CommonGenes_*.tsv", optional: true

  """
  Rscript ${scripts_dir}/differential_analysis/dexseq_sleuth_integration.R --genome_annotation ${genome_annotation} --sleuth_analysis ${sleuth_file} --dexseq_analysis ${dexseq_file} --p_thr ${params.pval_threshold}
  """

}