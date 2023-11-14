process AnalyzeSingleCellData {

  label 'python'

  publishDir "${projectDir}/${params.single_cell_analysis}", mode: "copy", pattern: "*.{pkl,png,h5,hdf5}"

  input:
  path scripts_dir
  path counts

  output:
  path "*.pkl"
  path "*.png"
  path "*.{h5,hdf5}"
  
  """
  python ${scripts_dir}/singel_cell/tenx_scrnaseq_analysis.py --min_counts 500 --min_detected_genes 500 --n_neighbors 5,30,50 --input_dir ./
  """

}