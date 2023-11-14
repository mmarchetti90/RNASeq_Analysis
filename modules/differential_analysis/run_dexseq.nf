process RunDEXSeq {

  label 'r'

  publishDir "${projectDir}/${params.dexseq_dir}", mode: 'copy', pattern: "DEU_*.{tsv,rds}"
  publishDir "${projectDir}/${params.dexseq_dir}/${params.dexseq_plots}", mode: 'copy', pattern: "DEU_*.{png}"

  input:
  each path(scripts_dir)
  each path(counts)
  path comparison

  output:
  path "DEU_*.rds", optional: true
  tuple val("${comparison.simpleName}"), path("DEU_*.tsv"), optional: false, emit: deu_files
  path "DEU_*_Dispersion-Plot.png", optional: true
  path "DEU_*_MA-Plot.png", optional: true

  """
  Rscript ${scripts_dir}/differential_analysis/run_dexseq.R --counts ${counts} --design ${comparison} --threads \$SLURM_CPUS_ON_NODE
  """

}