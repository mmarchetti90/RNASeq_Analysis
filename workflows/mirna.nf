
// ----------------Workflow---------------- //

include { ParseReadsList } from '../modules/parse_reads_list/parse_reads_list.nf'
include { GenerateBowtie2Index } from '../modules/bowtie/generate_bowtie_index.nf'
include { GetMirBaseAnnotation } from '../modules/mirbase/get_mirbase_annotation.nf'
include { TrimFastQ } from '../modules/trimgalore/trimgalore.nf'
include { RunBowtie2 } from '../modules/bowtie/run_bowtie.nf'
include { MergeCounts } from '../modules/merge_counts/merge_counts.nf'
include { RunDESeq2 } from '../modules/differential_analysis/run_deseq2.nf'

workflow MIRNA {

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
  
  // Download miRNA annotation from mirBase
  GetMirBaseAnnotation(genome_resources.flatten().first())

  // Creating channel for existing bowtie2 index, or building de novo
  if (new File("${params.bowtie_index_dir}/${params.bowtie_index_prefix}.1.bt2").exists()) {

    Channel
    .fromPath("${params.bowtie_index_dir}/*bt2")
    .collect()
    .set{ bowtie_index }

  }
  else {

    GenerateBowtie2Index(genome_fasta)
    bowtie_index = GenerateBowtie2Index.out.bowtie_index

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

  // BOWTIE2 ALIGNMENT -------------------- //

  // Run Bowtie2 alignment
  RunBowtie2(GenerateBowtie2Index.out.bowtie_index.collect().ifEmpty("$projectDir"), TrimFastQ.out.trimmed_fastq_files, GetMirBaseAnnotation.out.mirna_annotation)

  // Merge gene count files
  MergeCounts(RunBowtie2.out.mirna_counts.collect(), "miRNA")

  // DOWNSTREAM ANALISES ------------------ //

  if (params.alignment_only == false) {

    // DIFFERENTIAL EXPRESSION ANALYSIS --- //

    // Differential gene expression analysis with DESeq2
    RunDESeq2(scripts_dir, MergeCounts.out, comparison_files)

    // Collecting outputs as ready signal
    Channel.value("0")
      .join(MergeCounts.out.merged_counts.ifEmpty("0"))
      .join(RunDESeq2.out[0].ifEmpty("0"))
      .toList()
      .set{ done_mirna }

  }

}