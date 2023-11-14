process RunDESeq2 {

  label 'r'

  publishDir "${projectDir}/${params.deseq_dir}", mode: 'copy', pattern: "DEA_*.{tsv,rds}"
  publishDir "${projectDir}/${params.deseq_dir}/${params.deseq_plots}", mode: 'copy', pattern: "DEA_*.{png}"

  input:
  each path(scripts_dir)
  each path(counts)
  path comparison

  output:
  path "DEA_*.rds", optional: true
  path "DEA_*.tsv", optional: false, emit: dea_files
  path "DEA_*_MeanVsVariance.png", optional: true
  path "DEA_*_MA-Plot.png", optional: true
  path "DEA_*_SampleDistance.png", optional: true
  path "DEA_*_PCA.png", optional: true
  path "DEA_*_VolcanoPlot.png", optional: true

  """
  Rscript ${scripts_dir}/differential_analysis/run_deseq2.R --counts ${counts} --design ${comparison} --mincounts ${params.min_counts} --p_thr ${params.pval_threshold}
  """

}