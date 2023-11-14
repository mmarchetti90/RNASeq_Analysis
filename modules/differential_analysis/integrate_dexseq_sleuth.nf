process IntegrateAnalyses {

  label 'r'

  publishDir "${projectDir}/${params.dexseq_sleuth_integration}/${params.common_genes_dir}", mode: 'copy', pattern: "CommonGenes_*.tsv"
  publishDir "${projectDir}/${params.dexseq_sleuth_integration}/${params.enriched_terms_dir}", mode: 'copy', pattern: "Enrichment_*.tsv"
  publishDir "${projectDir}/${params.dexseq_sleuth_integration}/${params.enriched_go_terms_dir}", mode: 'copy', pattern: "Enrichment-GO*_*.tsv"
  publishDir "${projectDir}/${params.dexseq_sleuth_integration}/${params.clusterprofiler_rdata}", mode: 'copy', pattern: "clusterProfiler_*.RData.gz"

  input:
  each path(scripts_dir)
  each path(genome_annotation)
  tuple val(analysis_name), path(dexseq_file), path(sleuth_file)

  output:
  path "CommonGenes_*.tsv", optional: true
  path "Enrichment_*.tsv", optional: true
  path "Enrichment-GO*_*.tsv", optional: true
  path "clusterProfiler_*.RData.gz", optional: true

  """
  Rscript ${scripts_dir}/differential_analysis/dexseq_sleuth_integration.R --species ${params.species} --genome_annotation ${genome_annotation} --sleuth_analysis ${sleuth_file} --dexseq_analysis ${dexseq_file} --p_thr ${params.pval_threshold}
  """

}