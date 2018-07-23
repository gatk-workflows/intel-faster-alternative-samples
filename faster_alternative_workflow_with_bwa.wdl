workflow FasterAltWorkflow {

  File input1_fq
  File? input2_fq
  String sample
  String output_bam_basename
  String vcf_basename

  File wgs_calling_interval_list
  Int haplotype_scatter_count
  Int break_bands_at_multiples_of

  # Reference Files
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File ref_alt
  File ref_bwt
  File ref_sa
  File ref_amb
  File ref_ann
  File ref_pac

  # Path to tools
  String tool_path
  
  #Optimization flags
  Int compression_level
    
  call BamMapMarkSort {
   input: 
     input1_fq = input1_fq,
     input2_fq = input2_fq,
     sample = sample, 
     tool_path = tool_path,
     output_bam_basename = output_bam_basename,
     ref_dict = ref_dict,
     ref_fasta = ref_fasta,
     ref_fasta_index = ref_fasta_index,
     ref_amb = ref_amb,
     ref_ann = ref_ann,
     ref_bwt = ref_bwt,
     ref_pac = ref_pac,
     ref_sa = ref_sa
  }
  
  # Break the calling interval_list into sub-intervals
  # Perform variant calling on the sub-intervals, and then gather the results
  call ScatterIntervalList {
    input:
      interval_list = wgs_calling_interval_list,
      scatter_count = haplotype_scatter_count,
      break_bands_at_multiples_of = break_bands_at_multiples_of,
      compression_level = compression_level,
      tool_path = tool_path
  }
  
  # Call variants in parallel over WGS calling intervals
  scatter (index in range(ScatterIntervalList.interval_count)) {
    # Generate GVCF by interval
    call HaplotypeCaller {
      input:
        input_bam = BamMapMarkSort.output_bam,
        input_bam_index = BamMapMarkSort.output_bam_index,
        interval_list = ScatterIntervalList.out[index],
        vcf_basename = vcf_basename,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        compression_level = compression_level, 
        tool_path = tool_path
    }
  }

  # Combine by-interval GVCFs into a single sample GVCF file
  call GatherHC {
    input:
      input_vcfs = HaplotypeCaller.output_vcf,
      input_vcfs_indexes = HaplotypeCaller.output_vcf_index,
      vcf_basename = vcf_basename,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      tool_path = tool_path
  }
}

task BamMapMarkSort {
  File input1_fq
  File? input2_fq
  String sample
  String tool_path
  String bwa_options
  String output_bam_basename

  # reference files
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa

  # config
  Int cores
  String sort_memory
  Int cpu
  String memory
  String tmpDir
    
  command <<<
      set -euo pipefail
      ${tool_path}/bwa ${bwa_options} ${input1_fq} ${input2_fq} \
      | \
      ${tool_path}/samblaster \
      | \
      ${tool_path}/samtools sort -m ${sort_memory} -@ ${cores} -l 0 -o ${output_bam_basename}.bam /dev/stdin \
      && \
      ${tool_path}/samtools index ${output_bam_basename}.bam ${output_bam_basename}.bai
      #Using Samtools is 2x faster
      #${tool_path}/sambamba view -S -f bam -l 0 /dev/stdin \
      #| \
      #${tool_path}/sambamba sort --tmpdir=${tmpDir} -m ${sort_memory} -t ${cores} -l 0 -o ${output_bam_basename}.bam /dev/stdin
  >>>
  runtime {
    memory: memory
    cpu: cpu
  }
  output {
    File output_bam = "${output_bam_basename}.bam"
    File output_bam_index = "${output_bam_basename}.bai"
    #if using sambamba - use the .bam.bai index output file
    #File output_bam_index = "${output_bam_basename}.bam.bai"
  }
}

# This task calls picard's IntervalListTools to scatter the input interval list into scatter_count sub interval lists
# Note that the number of sub interval lists may not be exactly equal to scatter_count.  There may be slightly more or less.
# Thus we have the block of python to count the number of generated sub interval lists.
task ScatterIntervalList {
  File interval_list
  Int scatter_count
  Int break_bands_at_multiples_of
  Int compression_level
  String memory
  Int cpu
  String tool_path
  
  command <<<
    set -e
    mkdir out
    java -Dsamjdk.compression_level=${compression_level} -Xmx10g -jar ${tool_path}/picard.jar \
      IntervalListTools \
      SCATTER_COUNT=${scatter_count} \
      SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
      UNIQUE=true \
      SORT=true \
      BREAK_BANDS_AT_MULTIPLES_OF=${break_bands_at_multiples_of} \
      INPUT=${interval_list} \
      OUTPUT=out

    python3 <<CODE
    import glob, os
    # Works around a JES limitation where multiples files with the same name overwrite each other when globbed
    intervals = sorted(glob.glob("out/*/*.interval_list"))
    for i, interval in enumerate(intervals):
      (directory, filename) = os.path.split(interval)
      newName = os.path.join(directory, str(i + 1) + filename)
      os.rename(interval, newName)
    print(len(intervals))
    CODE
  >>>
  output {
    Array[File] out = glob("out/*/*.interval_list")
    Int interval_count = read_int(stdout())
  }
  runtime {
    memory: memory
    cpu: cpu
  }
}

# Call variants on a single sample with HaplotypeCaller to produce a GVCF
task HaplotypeCaller {
  File input_bam
  File input_bam_index
  File interval_list
  String vcf_basename
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  Float? contamination
  String gatk_gkl_pairhmm_implementation
  Int gatk_gkl_pairhmm_threads
  String smith_waterman_implementation
  Int compression_level
  String memory
  Int cpu
  String tool_path
  
  # We use interval_padding 500 below to make sure that the HaplotypeCaller has context on both sides around
  # the interval because the assembly uses them.
  command <<<
      export GATK_LOCAL_JAR=${tool_path}/gatk4.0.0.0/gatk.jar && \
      ${tool_path}/gatk4.0.0.0/gatk --java-options -Xmx6g \
      HaplotypeCaller \
      -R ${ref_fasta} \
      -I ${input_bam} \
      -O ${vcf_basename}.vcf.gz \
      -L ${interval_list} \
      -ip 100 \
      -contamination ${default=0 contamination} \
      --max-alternate-alleles 3 \
      -pairHMM ${gatk_gkl_pairhmm_implementation} \
      --native-pair-hmm-threads ${gatk_gkl_pairhmm_threads} \
      --smith-waterman ${smith_waterman_implementation}
      #-pairHMM AVX_LOGLESS_CACHING \
  >>>
  runtime {
    memory: memory
    cpu: cpu
    #require_fpga: "yes"
  }
  output {
    File output_vcf = "${vcf_basename}.vcf.gz"
    File output_vcf_index = "${vcf_basename}.vcf.gz.tbi"
  }
}

task GatherHC {
  Array[File] input_vcfs
  Array[File] input_vcfs_indexes
  String vcf_basename
  File ref_dict
  File ref_fasta
  File ref_fasta_index
  String tool_path
  String memory
  Int cpu

  command {
    java -cp ${tool_path}/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants \
      -V ${sep=' -V ' input_vcfs} \
      -R ${ref_fasta} \
      -out ${vcf_basename}.vcf \
      -assumeSorted
  }
  runtime {
    memory: memory
    cpu: cpu
  }
  output {
    File gathered_vcf = "${vcf_basename}.vcf"
    File gathered_vcf_index = "${vcf_basename}.vcf.idx"
  }
}
