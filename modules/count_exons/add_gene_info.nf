process AddGeneInfo {

    label 'python'

    publishDir "${projectDir}/${params.exon_counts_dir}", mode: "copy", pattern: "MergedExonCounts_Annotated.tsv"

    input:
    path scripts_dir
    path genome_annotation
    path counts_file

    output:
    path "MergedExonCounts_Annotated.tsv", emit: annotated_counts_file

    """
    python ${scripts_dir}/add_gene_info/exons_add_gene_info.py \
    --counts ${counts_file} \
    --gtf ${genome_annotation}

    mv annotated_counts.tsv MergedExonCounts_Annotated.tsv
    """

}