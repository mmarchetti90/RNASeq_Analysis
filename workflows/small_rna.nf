
// ----------------Workflow---------------- //

include { ParseReadsList } from '../modules/parse_reads_list/parse_reads_list.nf'
include { TrimFastQ } from '../modules/trimgalore/trimgalore.nf'
include { GenerateStarIndex } from '../modules/star/generate_star_index.nf'
include { RunSTAR } from '../modules/star/run_star_small_rna.nf'
include { MergeCounts } from '../modules/merge_counts/merge_counts.nf'
include { RunDESeq2 } from '../modules/differential_analysis/run_deseq2.nf'
include { FindEnrichment } from '../modules/cluster_profiler/find_enrichment.nf'

workflow SMALLRNA {

  main:
  // LOADING RESOURCES -------------------- //

  // Channel for the directory containing the scripts used by the pipeline
  Channel
    .fromPath("${projectDir}/scripts")
    .set{ scripts_dir }

  // Channel for genome fasta
  Channel
    .fromPath("${params.genome_fasta_path}")
    .set{ genome_fasta }

  // Channel for genome annotation
  Channel
    .fromPath("${params.genome_annotation_path}")
    .set{ genome_annotation }

  // Channel for transcripts fasta
  Channel
    .fromPath("${params.transcripts_fasta_path}")
    .set{ transcripts_fasta }

  // Creating channel for existing star index, or building de novo
  if (new File("${params.star_index_dir}/Genome").exists()) {

    Channel
    .fromPath("${params.star_index_dir}")
    .set{ star_index }

  }
  else {

    sjdboverhang = params.read_length - 1
    GenerateStarIndex(genome_fasta, genome_annotation, sjdboverhang)
    star_index = GenerateStarIndex.out.star_index

  }
  
  // Loading reads list file
  Channel
  .fromPath("${params.metadata_path}")
  .set{ metadata }

  // Parsing reads list to output lists fastq files, and downstream comparisons
  ParseReadsList(scripts_dir, metadata)

  // Creating raw_reads channel
  ParseReadsList.out.reads_list
    .splitCsv(header: true, sep: '\t')
    .map{row -> tuple(row.SampleID, file(row.File1), file(row.File2))}
    .set{ raw_reads }

  // Creating comparison files channel
  ParseReadsList.out.comparison_files
    .flatten()
    .set{ comparison_files }
  
  // TRIMGALORE --------------------------- //

  // Trimming adapters
  TrimFastQ(raw_reads)

  // STAR ALIGNMENT ----------------------- //

  // Run STAR alignment
  RunSTAR(star_index, TrimFastQ.out.trimmed_fastq_files)

  // Merge gene count files
  MergeCounts(RunSTAR.out.gene_counts.collect(), "Gene")

  // DOWNSTREAM ANALYSES ------------------ //

  if (params.alignment_only == false) {

    // DIFFERENTIAL EXPRESSION ANALYSIS --- //

    // Differential gene expression analysis with DESeq2
    RunDESeq2(scripts_dir, MergeCounts.out.merged_counts, comparison_files)

    // ENRICHMENT ANALYSIS ---------------- //

    // Enrichment analysis with clusterProfiler
    FindEnrichment(scripts_dir, RunDESeq2.out.dea_files)

  }

}