
// ----------------Workflow---------------- //

include { ParseReadsList } from '../modules/parse_reads_list/parse_reads_list.nf'
include { RunFastQC } from '../modules/fastqc/run_fastqc.nf'
include { RunCellranger } from '../modules/single_cell/run_cellranger.nf'
include { AnalyzeSingleCellData } from '../modules/single_cell/analyze_single_cell_data.nf'

workflow CELLRANGER {

  main:
  // LOADING RESOURCES -------------------- //

  // Channel for the directory containing the scripts used by the pipeline
  Channel
    .fromPath("${projectDir}/scripts")
    .set{ scripts_dir }

  // Channel for cellranger index
  Channel
    .fromPath("${params.cellranger_index_dir}")
    .set{ cellranger_index }

  // Loading reads list file
  Channel
  .fromPath("${params.metadata_path}")
  .set{ metadata }

  // Parsing reads list to output lists fastq files, and downstream comparisons
  ParseReadsList(scripts_dir, metadata)

  // Creating raw_reads channel
  ParseReadsList.out.reads_list
    .splitCsv(header: true, sep: '\t')
    .map{row -> tuple(row.SampleID, row.File1, file(row.File2))}
    .set{ raw_reads }

  // FASTQC ------------------------------- //

  // Process runs only if fastq file exists
  RunFastQC(raw_reads)

  // CELLRANGER RUN ----------------------- //

  // Cellranger alignment for 10X data
  RunCellranger(cellranger_index, raw_reads)

  // DOWNSTREAM ANALYSES ------------------ //

  if (params.alignment_only == false) {

    // SINGLE CELL ANALYSIS TESTING ----------------------- //

    // Generating UMAPs with 5, 30, and 50 neighbors. UMAP models are saved as Pickle files
    AnalyzeSingleCellData(scripts_dir, RunCellranger.out.cells_gene_counts.collect())

  }

}