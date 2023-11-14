process RunKallisto {

  label 'kallisto'

  publishDir "${projectDir}/${params.transcript_counts_dir}", mode: "copy", pattern: "*_kallisto_out"

  input:
  each path(index)
  tuple val(read_id), path(read1), path(read2)

  output:
  path "${read_id}_abundance.tsv", optional: true, emit: transcript_counts
  path "${read_id}_kallisto_out", optional: true, emit: kallisto_output
  
  """
  mkdir ${read_id}_kallisto_out

  if [[ "${read2}" == "mock.fastq" ]]
  then

    if [[ "${params.strand}" == "unstranded" ]]
    then

      kallisto quant \
      -i ${index} \
      -t \$SLURM_CPUS_ON_NODE \
      --single \
      -l 200 -s 20 \
      -b 100 \
      -o ${read_id}_kallisto_out \
      ${read1}

    elif [[ "${params.strand}" == "stranded" ]]
    then

      kallisto quant \
      -i ${index} \
      -t \$SLURM_CPUS_ON_NODE \
      --single \
      --fr-stranded \
      -l 200 -s 20 \
      -b 100 \
      -o ${read_id}_kallisto_out \
      ${read1}

    else

      kallisto quant \
      -i ${index} \
      -t \$SLURM_CPUS_ON_NODE \
      --single \
      --rf-stranded \
      -l 200 -s 20 \
      -b 100 \
      -o ${read_id}_kallisto_out \
      ${read1}

    fi

  else

    if [[ "${params.strand}" == "unstranded" ]]
    then

      kallisto quant \
      -i ${index} \
      -t \$SLURM_CPUS_ON_NODE \
      -b 100 \
      -o ${read_id}_kallisto_out \
      ${read1} ${read2}

    elif [[ "${params.strand}" == "stranded" ]]
    then

      kallisto quant \
      -i ${index} \
      -t \$SLURM_CPUS_ON_NODE \
      --fr-stranded \
      -b 100 \
      -o ${read_id}_kallisto_out \
      ${read1} ${read2}

    else

      kallisto quant \
      -i ${index} \
      -t \$SLURM_CPUS_ON_NODE \
      --rf-stranded \
      -b 100 \
      -o ${read_id}_kallisto_out \
      ${read1} ${read2}

    fi

  fi

  cp ${read_id}_kallisto_out/abundance.tsv ${read_id}_abundance.tsv
  """

}