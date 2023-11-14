
// ----------------Workflow---------------- //

include { ParseReadsList } from '../modules/parse_reads_list/parse_reads_list.nf'
include { RunFastQC } from '../modules/fastqc/run_fastqc.nf'
include { GenerateStarIndex } from '../modules/star/generate_star_index.nf'
include { RunSTAR } from '../modules/star/run_star_standard_rna.nf'
include { MergeCounts as MergeGeneCounts } from '../modules/merge_counts/merge_counts.nf'
include { RunDESeq2 } from '../modules/differential_analysis/run_deseq2.nf'
include { FindEnrichment } from '../modules/cluster_profiler/find_enrichment.nf'
include { CountExons } from '../modules/count_exons/count_exons.nf'
include { MergeCounts as MergeExonCounts } from '../modules/merge_counts/merge_counts.nf'
include { AddGeneInfo } from '../modules/count_exons/add_gene_info.nf'
include { RunDEXSeq } from '../modules/differential_analysis/run_dexseq.nf'
include { GenerateKallistoIndex } from '../modules/kallisto/generate_kallisto_index.nf'
include { RunKallisto } from '../modules/kallisto/run_kallisto.nf'
include { MergeCounts as MergeTranscriptCounts } from '../modules/merge_counts/merge_counts.nf'
include { RunSleuth } from '../modules/differential_analysis/run_sleuth.nf'
include { IntegrateAnalyses as IntegrateDEXSeqSleuth } from '../modules/differential_analysis/integrate_dexseq_sleuth.nf'

workflow STANDARD {

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

  // Creating channel for existing kallisto index, or building de novo
  if (new File("${params.kallisto_index_path}").exists()) {

    Channel
    .fromPath("${params.kallisto_index_path}")
    .set{ kallisto_index }

  }
  else {

    GenerateKallistoIndex(transcripts_fasta)
    kallisto_index = GenerateKallistoIndex.out.kallisto_index

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
  
  // FASTQC ------------------------------- //

  // Process runs only if fastq file exists
  RunFastQC(raw_reads)

  // STAR ALIGNMENT ----------------------- //

  // Run STAR alignment
  RunSTAR(star_index, raw_reads)

  // Merge gene count files
  MergeGeneCounts(RunSTAR.out.gene_counts.collect(), "Gene")

  // Counting exons
  CountExons(genome_annotation, RunSTAR.out.bam_files)

  // Merge exon count files
  MergeExonCounts(CountExons.out.exon_counts.collect(), "Exon")

  // Add gene info to exon counts
  AddGeneInfo(genome_annotation, MergeExonCounts.out.merged_counts)

  // KALLISTO ALIGNMENT ------------------- //

  // Run Kallisto pseudo-alignment
  RunKallisto(kallisto_index, raw_reads)

  // Merge trascript count files
  MergeTranscriptCounts(RunKallisto.out.transcript_counts.collect(), "Transcript")

  // DOWNSTREAM ANALYSES ------------------ //

  if (params.alignment_only == false) {

    // DIFFERENTIAL GENE EXPRESSION ------- //

    // Differential gene expression analysis with DESeq2
    RunDESeq2(scripts_dir, MergeGeneCounts.out.merged_counts, comparison_files)

    // GENE LEVEL ENRICHMENT ANALYSIS ----- //

    // Enrichment analysis with clusterProfiler
    FindEnrichment(scripts_dir, RunDESeq2.out.dea_files)

    // DIFFERENTIAL EXON USAGE ------------ //

    // Differential exon utilization analysis with DEXSeq
    RunDEXSeq(scripts_dir, AddGeneInfo.out.annotated_counts_file, comparison_files)

    // DIFFERENTIAL TRANSCRIPT EXPRESSION - //

    // Differential transcript expression analysis with Sleuth
    RunSleuth(scripts_dir, RunKallisto.out.kallisto_output.collect(), comparison_files)

    // DEXSEQ SLEUTH INTEGRATION ---------- //

    // Merging necessary channels
    RunDEXSeq.out.deu_files
      .join(RunSleuth.out.dea_files, by: 0, remainder: false)
      .set{ analysis_integration_input }

    // Integrating DEXSeq and Sleuth results, running clusterProfiler on genes with both differentially used exons and differentially expressed transcripts
    IntegrateDEXSeqSleuth(scripts_dir, genome_annotation, analysis_integration_input)

  }

}