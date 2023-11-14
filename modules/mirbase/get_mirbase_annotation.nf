process GetMirBaseAnnotation {

  label 'python'

  output:
  path "*gff3", optional: true, emit: mirna_annotation

  """
  if [ "${params.species}" == "Caenorhabditis_elegans" ]
  then

    curl --output cel.gff3 https://www.mirbase.org/ftp/CURRENT/genomes/cel.gff3

  elif [ "${params.species}" == "Drosophila_melanogaster" ]
  then

    curl --output dme.gff3 https://www.mirbase.org/ftp/CURRENT/genomes/dme.gff3

  elif [ "${params.species}" == "Danio_rerio" ]
  then

    curl --output dre.gff3 https://www.mirbase.org/ftp/CURRENT/genomes/dre.gff3

  elif [ "${params.species}" == "Homo_sapiens" ]
  then

    curl --output hsa.gff3 https://www.mirbase.org/ftp/CURRENT/genomes/hsa.gff3

  elif [ "${params.species}" == "Mus_musculus" ]
  then

    curl --output mmu.gff3 https://www.mirbase.org/ftp/CURRENT/genomes/mmu.gff3

  elif [ "${params.species}" == "Pan_troglodytes" ]
  then

    curl --output ptr.gff3 https://www.mirbase.org/ftp/CURRENT/genomes/ptr.gff3

  elif [ "${params.species}" == "Rattus_norvegicus" ]
  then

    curl --output rno.gff3 https://www.mirbase.org/ftp/CURRENT/genomes/rno.gff3

  elif [ "${params.species}" == "Xenopus_tropicalis" ]
  then

    curl --output xtr.gff3 https://www.mirbase.org/ftp/CURRENT/genomes/xtr.gff3

  else

    touch mock.txt

  fi
  """

}