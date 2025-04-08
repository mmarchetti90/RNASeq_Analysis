#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
Pipeline to analyze bulk and single-cell mRNA-Seq data for gene- and transcript-level analyses
*/

// ----------------Workflow---------------- //

include { STANDARD } from './workflows/standard.nf'
include { GENELVL } from './workflows/gene_level.nf'
include { TRANSCRIPTLVL } from './workflows/transcript_level.nf'
include { CELLRANGER } from './workflows/cellranger.nf'
include { SMALLRNA } from './workflows/small_rna.nf'
include { MIRNA } from './workflows/mirna.nf'

workflow {

  // WORKFLOW SELECTION ------------------- //

  if ("$workflow.profile" == "standard") {
    
    STANDARD()

  }
  else if ("$workflow.profile" == "gene_lvl") {
    
    GENELVL()

  }
  else if ("$workflow.profile" == "transcript_lvl") {
    
    TRANSCRIPTLVL()

  }
  else if ("$workflow.profile" == "cellranger") {
    
    CELLRANGER()

  }
  else if ("$workflow.profile" == "small_rna") {
    
    SMALLRNA()

  }
  else if ("$workflow.profile" == "mirna") {
    
    MIRNA()

  } else {

    println "ERROR: Unrecognized profile!"
    println "Please chose one of: standard, gene_lvl, transcript_lvl, cellranger, small_rna, mirna"

  }

}