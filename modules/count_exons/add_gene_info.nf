process AddGeneInfo {

    label 'subread'

    publishDir "${projectDir}/${params.exon_counts_dir}", mode: "copy", pattern: "MergedExonCounts_Annotated.tsv"

    input:
    path genome_annotation
    path counts_file

    output:
    path "MergedExonCounts_Annotated.tsv", emit: annotated_counts_file

    """
    # Making ExonID to GeneID conversion table
    echo "Making ExonID to GeneID conversion table"

    awk '\$3 == "exon" { print }' ${genome_annotation} | \
    grep -o -e 'gene_id ".*"' | \
    sed 's/.........//' | \
    sed 's/".*//g' > gene_list.txt

    awk '\$3 == "exon" { print }' ${genome_annotation} | \
    grep -o -e 'exon_id ".*"' | \
    sed 's/.........//' | \
    sed 's/".*//g' > exon_list.txt

    echo -e "ExonID\tGeneID" > conversion_table.txt
    paste exon_list.txt gene_list.txt | sort | uniq >> conversion_table.txt

    rm gene_list.txt exon_list.txt

    # Merging to counts file
    echo "Adding GeneID info"
    paste conversion_table.txt ${counts_file} | cut -f1-2,4- > MergedExonCounts_Annotated.tsv

    rm conversion_table.txt
    """

}