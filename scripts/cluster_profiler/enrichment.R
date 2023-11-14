#!/usr/bin/Rscript

# This script reads in differentially expressed genes (DEG) and performs enrichment analyses

### ---------------------------------------- ###

parseArgs <- function() {
  
  # Read command line arguments
  args <- commandArgs()
  
  # DEA file location
  filepath <- args[match("--filepath", args) + 1]
  
  # Species
  species <- args[match("--species", args) + 1]
  
  # Log2FC threshold
  if("--log2fc_thr" %in% args) {
    
    log2fc <- args[match("--log2fc_thr", args) + 1]
    
  } else {
    
    log2fc <- 0
    
  }
  
  # Adjusted p value threshold
  if("--p_thr" %in% args) {
    
    pval <- args[match("--p_thr", args) + 1]
    
  } else {
    
    pval <- 0.05
    
  }
  
  # msigdbr categories file
  if("--msigdbr_categories_file" %in% args) {
    
    msigdbr_categories_file <- args[match("--msigdbr_categories_file", args) + 1]
    
  } else {
    
    msigdbr_categories_file <- ""
    
  }
  
  return(c(filepath, species, log2fc, pval, msigdbr_categories_file))
  
}

### ---------------------------------------- ###

loadOrgDB <- function(org) {
  
  if(org == "Homo_sapiens") {
    library(org.Hs.eg.db)
    return(org.Hs.eg.db)
  } else if(org == "Mus_musculus") {
    library(org.Mm.eg.db)
    return(org.Mm.eg.db)
  } else if(org == "Rattus_norvegicus") {
    library(org.Rn.eg.db)
    return(org.Rn.eg.db)
  } else if(org == "Danio_rerio") {
    library(org.Dr.eg.db)
    return(org.Dr.eg.db)
  } else if(org == "Drosophila_melanogaster") {
    library(org.Dm.eg.db)
    return(org.Dm.eg.db)
  } else if(org == "Saccharomyces_cerevisiae") {
    library(org.Sc.eg.db)
    return(org.Sc.eg.db)
  } else if(org == "Caenorhabditis_elegans") {
    library(org.Ce.eg.db)
    return(org.Ce.eg.db)
  } else  {
    print("Unrecognized species...")
    quit(save = "no")
  }
  
}

### ---------------------------------------- ###

importDEA <- function(params, ags) {
  
  dea <- read.delim(as.character(params[1]), header = TRUE, row.names = 1, sep = "\t")
  
  if(nrow(dea) == 0) {
    
    print("DEA file is empty...")
    quit(save = "no")
    
  }
  
  # Converting gene names to Ensemble (if not Ensembl already*)
  # (* check only if the top 10 genes have an Ensembl name)
  symbol_score <- sum(tolower(rownames(dea)[1:10]) %in% tolower(ags$gene_symbol))
  entrez_score <- sum(rownames(dea)[1:10] %in% ags$entrez_gene)
  ensemble_score <- sum(rownames(dea)[1:10] %in% ags$ensembl_gene)
  
  if(entrez_score > ensemble_score & entrez_score > symbol_score) {
    
    conversion <- unique(ags$ensembl_gene)[match(rownames(dea), unique(ags$entrez_gene))]
    unknown <- is.na(conversion)
    conversion[unknown] <- rownames(dea)[unknown]
    conversion[duplicated(conversion)] <- rownames(dea)[duplicated(conversion)]
    rownames(dea) <- conversion
    
  }
  
  if(symbol_score > ensemble_score & symbol_score > entrez_score) {
    
    conversion <- unique(ags$ensembl_gene)[match(tolower(rownames(dea)), tolower(unique(ags$gene_symbol)))]
    unknown <- is.na(conversion)
    conversion[unknown] <- rownames(dea)[unknown]
    conversion[duplicated(conversion)] <- rownames(dea)[duplicated(conversion)]
    rownames(dea) <- conversion
    
  }
  
  # Subset DEA into AllGenes, UpregGenes, and DownregGenes, while filtering for log2fc and padj (or pvalue if less than 10 genes have padj < p_thr)
  if(NROW(subset(dea, abs(dea$log2FoldChange) >= as.numeric(params[3]) & dea$padj < as.numeric(params[4]))) >= 10) {

    datasets <- list("AllGenes" = subset(dea, abs(dea$log2FoldChange) >= as.numeric(params[3]) & dea$padj < as.numeric(params[4])),
                     "UpregGenes" = subset(dea, dea$log2FoldChange >= as.numeric(params[3]) & dea$padj < as.numeric(params[4])),
                     "DownregGenes" = subset(dea, dea$log2FoldChange <= -as.numeric(params[3]) & dea$padj < as.numeric(params[4])))

  } else {

    datasets <- list("AllGenes" = subset(dea, abs(dea$log2FoldChange) >= as.numeric(params[3]) & dea$pvalue < as.numeric(params[4])),
                     "UpregGenes" = subset(dea, dea$log2FoldChange >= as.numeric(params[3]) & dea$pvalue < as.numeric(params[4])),
                     "DownregGenes" = subset(dea, dea$log2FoldChange <= -as.numeric(params[3]) & dea$pvalue < as.numeric(params[4])))

  }
  
  return(datasets)
  
}

### ---------------------------------------- ###

runEnricher <- function(out, ds, t2g) {
  
  # Run analysis
  enrichment <- enricher(gene = rownames(ds), TERM2GENE = t2g, pvalueCutoff = 0.1, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 1000)
  
  if(! is.null(enrichment)) {
    
    # Export file
    write.table(enrichment, out, row.names = FALSE, sep='\t')
    
  }
  
  return(enrichment)
  
}

### ---------------------------------------- ###

runEnrichGO <- function(out, odb, ont, ds) {
  
  # Run analysis
  enrichment <- enrichGO(gene = rownames(ds), ont = ont, OrgDb = odb, keyType = "ENSEMBL", pvalueCutoff = 0.1, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 1000)
  
  if(! is.null(enrichment)) {
    
    # Simplifying
    enrichment <- clusterProfiler::simplify(enrichment, cutoff = 0.7)
    
    # Export file
    write.table(enrichment, out, row.names = FALSE, sep='\t')
    
  }
  
  return(enrichment)
  
}

### ---------------------------------------- ###

runGSEGO <- function(out, odb, ont, ds) {
  
  # Creating gene list
  gene_list <- ds$log2FoldChange
  names(gene_list) <- as.character(rownames(ds))
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # Run analysis
  enrichment <- gseGO(gene_list, ont = ont, OrgDb = odb, keyType = "ENSEMBL", pvalueCutoff = 0.1, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 1000, by = "fgsea")
  
  if(nrow(enrichment) != 0) {
    
    # Simplifying
    enrichment <- clusterProfiler::simplify(enrichment, cutoff = 0.7)
    
    # Export file
    write.table(enrichment, out, row.names = FALSE, sep='\t')
    
  }
  
  # Plot top 3 terms
  nterms <- min(3, nrow(enrichment))
  if(nterms > 0) {
    
    out <- paste(substring(out, 0, nchar(out) - 4), ".png", sep = "")
    png(file = out, width = 1536, height = 512)
    gsea_plot <- gseaplot2(enrichment, geneSetID = enrichment$ID[1 : nterms], pvalue_table = T, base_size = 15)
    print(gsea_plot)
    dev.off()
    
  }
  
  return(enrichment)
  
}

### ---------------------------------------- ###

runGSEA <- function(out, ds, t2g) {
  
  # Creating gene list
  gene_list <- ds$log2FoldChange
  names(gene_list) <- as.character(rownames(ds))
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # Run analysis
  enrichment <- GSEA(gene_list, pvalueCutoff = 0.05, pAdjustMethod = "BH", eps = 0, by = "fgsea", TERM2GENE = t2g)
  
  if(nrow(enrichment) != 0) {
    
    # Export file
    write.table(enrichment, out, row.names = FALSE, sep='\t')
    
  }
  
  # Plot top 3 terms
  nterms <- min(3, nrow(enrichment))
  if(nterms > 0) {
    
    out <- paste(substring(out, 0, nchar(out) - 4), ".png", sep = "")
    png(file = out, width = 1536, height = 512)
    gsea_plot <- gseaplot2(enrichment, geneSetID = enrichment$ID[1 : nterms], pvalue_table = T, base_size = 15)
    print(gsea_plot)
    dev.off()
    
  }
  
  return(enrichment)
  
}

### ------------------MAIN------------------ ###

library(msigdbr)
library(clusterProfiler)
library(enrichplot)

parameters <- parseArgs()
species <- paste(strsplit(parameters[2], "_")[[1]], collapse = " ")

# Loading organism database
org_db <- loadOrgDB(parameters[2])

# Loading gene sets (all of them)
if(species %in% msigdbr_species()$species_name) {
  
  all_gene_sets <- msigdbr(species = species)
  
} else {
  
  print("Unrecognized species...")
  quit(save = "no")
  
}

# Filtering gene sets for user-specified ones
if(parameters[5] != "") {
  
  msigdbr_categories <- as.vector(read.table(parameters[5], sep='\t')[,1])
  all_gene_sets <- all_gene_sets[(all_gene_sets$gs_cat %in% msigdbr_categories) |
                                 (all_gene_sets$gs_subcat %in% msigdbr_categories) |
                                 (all_gene_sets$gs_name %in% msigdbr_categories),]
  
}

# Import DEA and subset into AllGenes, UpregGenes, and DownregGenes, while filtering for log2fc and padj
datasets <- importDEA(parameters, all_gene_sets)

# Creating term2gene matrix for enricher, but removing GO terms. They'll be processed and simplified separately
term2gene <- all_gene_sets[!(all_gene_sets$gs_subcat %in% c("GO:BP", "GO:MF", "GO:CC")), c("gs_name", "ensembl_gene")]

# Init lists of enrichment/gsea results to be saved as an RData object
enrichment_analyses <- list()
gsea_analyses <- list()

# Enrichment analysis
for(name in names(datasets)) {
  
  if(nrow(datasets[[name]]) == 0) {
    
    next
    
  }
  
  # Enricher
  output_name <- paste("Enrichment_", name, "_", gsub("DEA_", "", parameters[1]), sep = "")
  enrichment_data <- runEnricher(output_name, datasets[[name]], term2gene)
  
  data_name <- paste("Enrichment_AllSets_", name, sep = "")
  enrichment_analyses[[data_name]] <- enrichment_data
  
  # EnrichGO
  for(ontology in c("BP", "MF", "CC")) {
    
    output_name <- paste("Enrichment-GO", ontology, "_", name, "_", gsub("DEA_", "", parameters[1]), sep = "")
    runEnrichGO(output_name, org_db, ontology, datasets[[name]])
    
    data_name <- paste("Enrichment_", ontology, "_", name, sep = "")
    enrichment_analyses[[data_name]] <- enrichment_data
    
  }
  
}

# GO-term GSEA
for(ontology in c("BP", "MF", "CC")) {
  
  if(nrow(datasets$AllGenes) == 0) {
    
    next
    
  }
  
  output_name <- paste("GSEA-GO", ontology, "_AllGenes_", gsub("DEA_", "", parameters[1]), sep = "")
  gsea_data <- runGSEGO(output_name, org_db, ontology, datasets$AllGenes)
  
  if(nrow(gsea_data) != 0) {
    
    data_name <- paste("GSEA-GO", ontology, "_AllGenes", sep = "")
    gsea_analyses[[data_name]] <- gsea_data
    
  }
  
}

# GSEA analysis
if(nrow(datasets$AllGenes) != 0) {
  
  output_name <- paste("GSEA_AllGenes_", gsub("DEA_", "", parameters[1]), sep = "")
  gsea_data <- runGSEA(output_name, datasets$AllGenes, term2gene)
  
  if(nrow(gsea_data) != 0) {
    
    data_name <- "GSEA_AllSets_AllGenes"
    gsea_analyses[[data_name]] <- gsea_data
    
  }
  
}

# Saving datasets to RData object
sample_name <- gsub("DEA_", "", parameters[1])
output_name <- paste("clusterProfiler_", substring(sample_name, 1, nchar(sample_name) - 4), '.RData.gz', sep = "")
save(parameters, species, all_gene_sets, datasets, term2gene, enrichment_analyses, gsea_analyses, file = output_name, compress = T)
