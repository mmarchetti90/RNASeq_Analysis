process GenerateStarIndex {

  label 'star'

  publishDir "${projectDir}/${params.resources_dir}", mode: "copy", pattern: "${params.star_index_subdir}"

  input:
  path genome_fasta
  path genome_annotation
  val overhang

  output:
  path "${params.star_index_subdir}", emit: star_index

  script:
  """
  mkdir ${params.star_index_subdir}

  STAR \
  --runThreadN \$SLURM_CPUS_ON_NODE \
  --runMode genomeGenerate \
  --genomeDir ./${params.star_index_subdir}/ \
  --genomeFastaFiles ${genome_fasta} \
  --sjdbGTFfile ${genome_annotation} \
  --sjdbOverhang ${overhang} \
  --genomeSAindexNbases ${params.saindexnbases}
  """

}