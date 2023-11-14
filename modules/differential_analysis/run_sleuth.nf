process RunSleuth {

  label 'r'

  publishDir "${projectDir}/${params.sleuth_dir}", mode: 'copy', pattern: "Sleuth_*.{tsv,rds}"
  publishDir "${projectDir}/${params.sleuth_dir}/${params.sleuth_plots}", mode: 'copy', pattern: "Sleuth_*.{png}"

  input:
  each path(scripts_dir)
  path kallisto_dirs
  path comparison

  output:
  path "Sleuth_*.rds", optional: true
  tuple val("${comparison.simpleName}"), path("Sleuth_LRT_*.tsv"), optional: false, emit: dea_files
  path "Sleuth_Wald_*.tsv", optional: true
  path "Sleuth_*.png", optional: true

  """
  Rscript ${scripts_dir}/differential_analysis/run_sleuth.R --kallisto_out_path . --design ${comparison} --num_cores \$SLURM_CPUS_ON_NODE
  """

}