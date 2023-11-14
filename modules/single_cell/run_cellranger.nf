process RunCellranger {
	
  label 'cellranger'

  publishDir "${projectDir}/${params.cellranger_ouput_dir}", mode: "copy", pattern: "${file_id}"
  publishDir "${projectDir}/${params.cellranger_ouput_dir}", mode: "copy", pattern: "*_filtered_feature_bc_matrix.h5"

  input:
  each path(cellranger_index)
  tuple val(file_id), val(sample_prefix), path(sample_dir)

  output:
  path "${file_id}_filtered_feature_bc_matrix.h5", optional: true, emit: cells_gene_counts
  path "${file_id}", optional: true, emit: cellranger_output

  """
  cellranger count \
  --id=${file_id} \
  --transcriptome=${cellranger_index} \
  --fastqs=${sample_dir} \
  --sample=${sample_prefix} \
  --chemistry=${params.chemistry} \
  --nosecondary \
  --no-bam

  # Renaming the counts file
  mv ./${file_id}/outs/filtered_feature_bc_matrix.h5 ./${file_id}_filtered_feature_bc_matrix.h5
  """

}