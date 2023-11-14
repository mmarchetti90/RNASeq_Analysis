process FindEnrichment {

  label 'r'

  publishDir "${projectDir}/${params.enrichment_dir}/${params.enriched_terms_dir}", mode: 'copy', pattern: "Enrichment_*.tsv"
  publishDir "${projectDir}/${params.enrichment_dir}/${params.enriched_go_terms_dir}", mode: 'copy', pattern: "Enrichment-GO*_*.tsv"
  publishDir "${projectDir}/${params.enrichment_dir}/${params.gsea_terms_dir}", mode: 'copy', pattern: "GSEA_AllGenes_*.tsv"
  publishDir "${projectDir}/${params.enrichment_dir}/${params.gsea_terms_dir}/${params.gsea_terms_plots}", mode: 'copy', pattern: "GSEA_AllGenes_*.png"
  publishDir "${projectDir}/${params.enrichment_dir}/${params.gsea_go_terms_dir}", mode: 'copy', pattern: "GSEA-GO*_AllGenes_*.tsv"
  publishDir "${projectDir}/${params.enrichment_dir}/${params.gsea_go_terms_dir}/${params.gsea_go_terms_plots}", mode: 'copy', pattern: "GSEA-GO*_AllGenes_*.png"
  publishDir "${projectDir}/${params.enrichment_dir}/${params.clusterprofiler_rdata}", mode: 'copy', pattern: "clusterProfiler_*.RData.gz"

  input:
  each path(scripts_dir)
  path dea

  output:
  path "Enrichment_*.tsv", optional: true
  path "Enrichment-GO*_*.tsv", optional: true
  path "GSEA_AllGenes_*.tsv", optional: true
  path "GSEA_AllGenes_*.png", optional: true
  path "GSEA-GO*_AllGenes_*.tsv", optional: true
  path "GSEA-GO*_AllGenes_*.png", optional: true
  path "clusterProfiler_*.RData.gz", optional: true

  """
  Rscript ${scripts_dir}/cluster_profiler/enrichment.R --filepath ${dea} --species ${params.species} --log2fc_thr ${params.log2fc_threshold} --p_thr ${params.pval_threshold}
  """

}