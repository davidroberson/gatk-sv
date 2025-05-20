#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: Workflow
label: GatherSampleEvidence
doc: |
  GATK-SV workflow to gather evidence for structural variants
  from a single BAM/CRAM file. This runs various SV callers
  including Manta, MELT, Scramble, and Whamg.

requirements:
  SchemaDefRequirement:
    types:
      - $import: Structs.yml
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
  SubworkflowFeatureRequirement: {}
  MultipleInputFeatureRequirement: {}

inputs:
  # Required inputs
  bam_or_cram_file:
    type: File
    secondaryFiles: [.bai, .crai]
    doc: "Input BAM/CRAM file"
  
  bam_or_cram_index:
    type: File?
    doc: "Index for BAM/CRAM file (optional if using standard naming)"
  
  sample_id:
    type: string
    doc: "Sample ID"
  
  is_dragen_3_7_8:
    type: boolean?
    doc: "Whether alignment was done with Dragen 3.7.8"
  
  # Evidence collection flags
  collect_coverage:
    type: boolean
    default: true
    doc: "Whether to collect coverage data"
  
  collect_pesr:
    type: boolean
    default: true
    doc: "Whether to collect paired-end and split-read data"
  
  move_bam_or_cram_files:
    type: boolean
    default: false
    doc: "Whether to move BAM/CRAM files instead of copying them"
  
  run_localize_reads:
    type: boolean
    default: true
    doc: "Whether to localize reads first"
  
  # Common parameters
  primary_contigs_list:
    type: File
    doc: "List of primary contigs"
  
  reference_fasta:
    type: File
    secondaryFiles: [.fai, ^.dict]
    doc: "Reference genome FASTA"
  
  reference_index:
    type: File
    doc: "Reference genome FASTA index (.fai)"
  
  reference_dict:
    type: File
    doc: "Reference genome dictionary (.dict)"
  
  reference_version:
    type: string?
    doc: "Reference genome version ('38' or '19')"
  
  # BWA reference index files (for Dragen 3.7.8)
  reference_bwa_alt:
    type: File?
    doc: "BWA ALT file"
  
  reference_bwa_amb:
    type: File?
    doc: "BWA AMB file"
  
  reference_bwa_ann:
    type: File?
    doc: "BWA ANN file"
  
  reference_bwa_bwt:
    type: File?
    doc: "BWA BWT file"
  
  reference_bwa_pac:
    type: File?
    doc: "BWA PAC file"
  
  reference_bwa_sa:
    type: File?
    doc: "BWA SA file"
  
  # Coverage collection inputs
  preprocessed_intervals:
    type: File
    doc: "Preprocessed intervals for coverage collection"
  
  mem_gb_for_collect_counts:
    type: float?
    doc: "Memory in GB for collect counts task"
  
  disk_space_gb_for_collect_counts:
    type: int?
    doc: "Disk space in GB for collect counts task"
  
  # Manta inputs
  manta_region_bed:
    type: File
    doc: "BED file of regions for Manta"
  
  manta_region_bed_index:
    type: File
    doc: "Index for BED file of regions for Manta"
  
  manta_jobs_per_cpu:
    type: float?
    doc: "Manta jobs per CPU"
  
  manta_mem_gb_per_job:
    type: int?
    doc: "Manta memory in GB per job"
  
  # PESR inputs
  sd_locs_vcf:
    type: File
    doc: "VCF of segmental duplication locations"
  
  # Melt inputs
  melt_standard_vcf_header:
    type: File?
    doc: "Standard VCF header for MELT"
  
  melt_metrics_intervals:
    type: File?
    doc: "Intervals for MELT metrics"
  
  insert_size:
    type: float?
    doc: "Insert size for MELT"
  
  read_length:
    type: int?
    doc: "Read length for MELT"
  
  coverage:
    type: float?
    doc: "Coverage for MELT"
  
  metrics_intervals:
    type: File?
    doc: "Intervals for metrics collection"
  
  pct_chimeras:
    type: float?
    doc: "Percent chimeras for MELT"
  
  total_reads:
    type: float?
    doc: "Total reads for MELT"
  
  pf_reads_improper_pairs:
    type: int?
    doc: "PF reads improper pairs for MELT"
  
  # Scramble inputs
  mei_bed:
    type: File
    doc: "BED file of mobile element insertions for Scramble"
  
  scramble_alignment_score_cutoff:
    type: int?
    doc: "Alignment score cutoff for Scramble"
  
  scramble_percent_align_cutoff:
    type: int?
    doc: "Percent alignment cutoff for Scramble"
  
  scramble_min_clipped_reads_fraction:
    type: float?
    doc: "Minimum clipped reads fraction for Scramble"
  
  scramble_part2_threads:
    type: int?
    doc: "Threads for Scramble part 2"
  
  scramble_vcf_script:
    type: File?
    doc: "Script for Scramble VCF creation"
  
  # Required if running Scramble but not running Manta
  manta_vcf_input:
    type: File?
    doc: "Input Manta VCF (if not running Manta)"
  
  manta_vcf_index_input:
    type: File?
    doc: "Input Manta VCF index (if not running Manta)"
  
  # Wham inputs
  wham_include_list_bed_file:
    type: File
    doc: "BED file of regions to include for Wham"
  
  # Module metrics parameters
  run_module_metrics:
    type: boolean
    default: true
    doc: "Whether to run module metrics"
  
  primary_contigs_fai:
    type: File?
    doc: "FAI file for primary contigs (required if run_module_metrics is true)"
  
  baseline_manta_vcf:
    type: File?
    doc: "Baseline Manta VCF for metrics"
  
  baseline_wham_vcf:
    type: File?
    doc: "Baseline Wham VCF for metrics"
  
  baseline_melt_vcf:
    type: File?
    doc: "Baseline MELT VCF for metrics"
  
  baseline_scramble_vcf:
    type: File?
    doc: "Baseline Scramble VCF for metrics"
  
  # Docker
  sv_pipeline_docker:
    type: string
    doc: "SV pipeline Docker image"
  
  sv_base_mini_docker:
    type: string
    doc: "SV base mini Docker image"
  
  samtools_cloud_docker:
    type: string
    doc: "Samtools cloud Docker image"
  
  manta_docker:
    type: string?
    doc: "Manta Docker image"
  
  melt_docker:
    type: string?
    doc: "MELT Docker image"
  
  scramble_docker:
    type: string?
    doc: "Scramble Docker image"
  
  wham_docker:
    type: string?
    doc: "Wham Docker image"
  
  gatk_docker:
    type: string
    doc: "GATK Docker image"
  
  gatk_docker_pesr_override:
    type: string?
    doc: "GATK Docker image override for PESR"
  
  genomes_in_the_cloud_docker:
    type: string
    doc: "Genomes in the cloud Docker image"
  
  cloud_sdk_docker:
    type: string
    doc: "Cloud SDK Docker image"
  
  # Runtime attributes
  runtime_attr_localize_reads:
    type: Structs.yml#RuntimeAttr?
    doc: "Runtime attributes for localize reads"
  
  runtime_attr_split_cram:
    type: Structs.yml#RuntimeAttr?
    doc: "Runtime attributes for split CRAM"
  
  runtime_attr_concat_bam:
    type: Structs.yml#RuntimeAttr?
    doc: "Runtime attributes for concatenate BAM"
  
  runtime_attr_manta:
    type: Structs.yml#RuntimeAttr?
    doc: "Runtime attributes for Manta"
  
  runtime_attr_melt_coverage:
    type: Structs.yml#RuntimeAttr?
    doc: "Runtime attributes for MELT coverage"
  
  runtime_attr_melt_metrics:
    type: Structs.yml#RuntimeAttr?
    doc: "Runtime attributes for MELT metrics"
  
  runtime_attr_melt:
    type: Structs.yml#RuntimeAttr?
    doc: "Runtime attributes for MELT"
  
  runtime_attr_scramble_part1:
    type: Structs.yml#RuntimeAttr?
    doc: "Runtime attributes for Scramble part 1"
  
  runtime_attr_scramble_part2:
    type: Structs.yml#RuntimeAttr?
    doc: "Runtime attributes for Scramble part 2"
  
  runtime_attr_scramble_make_vcf:
    type: Structs.yml#RuntimeAttr?
    doc: "Runtime attributes for Scramble make VCF"
  
  runtime_attr_realign_soft_clips:
    type: Structs.yml#RuntimeAttr?
    doc: "Runtime attributes for realign soft clips"
  
  runtime_attr_scramble_part1_realigned:
    type: Structs.yml#RuntimeAttr?
    doc: "Runtime attributes for Scramble part 1 realigned"
  
  runtime_attr_scramble_part2_realigned:
    type: Structs.yml#RuntimeAttr?
    doc: "Runtime attributes for Scramble part 2 realigned"
  
  runtime_attr_scramble_make_vcf_realigned:
    type: Structs.yml#RuntimeAttr?
    doc: "Runtime attributes for Scramble make VCF realigned"
  
  runtime_attr_pesr:
    type: Structs.yml#RuntimeAttr?
    doc: "Runtime attributes for PESR"
  
  runtime_attr_wham:
    type: Structs.yml#RuntimeAttr?
    doc: "Runtime attributes for Wham"

# TODO: Complete the actual workflow steps

outputs:
  coverage_counts:
    type: File?
    doc: "Coverage counts file"
    # outputSource: CollectCounts/counts
  
  manta_vcf:
    type: File?
    secondaryFiles: [.tbi]
    doc: "Manta VCF"
    # outputSource: [Manta/vcf, manta_vcf_input]
  
  manta_index:
    type: File?
    doc: "Manta VCF index"
    # outputSource: [Manta/index, manta_vcf_index_input]
  
  melt_vcf:
    type: File?
    secondaryFiles: [.tbi]
    doc: "MELT VCF"
    # outputSource: MELT/vcf
  
  melt_index:
    type: File?
    doc: "MELT VCF index"
    # outputSource: MELT/index
  
  melt_coverage:
    type: float?
    doc: "MELT coverage"
    # outputSource: MELT/coverage_out
  
  melt_read_length:
    type: int?
    doc: "MELT read length"
    # outputSource: MELT/read_length_out
  
  melt_insert_size:
    type: float?
    doc: "MELT insert size"
    # outputSource: MELT/insert_size_out
  
  scramble_vcf:
    type: File?
    secondaryFiles: [.tbi]
    doc: "Scramble VCF"
    # outputSource: [ScrambleRealigned/vcf, Scramble/vcf]
  
  scramble_index:
    type: File?
    doc: "Scramble VCF index"
    # outputSource: [ScrambleRealigned/index, Scramble/index]
  
  scramble_clusters:
    type: File?
    doc: "Scramble clusters"
    # outputSource: [ScrambleRealigned/clusters, Scramble/clusters]
  
  scramble_table:
    type: File?
    doc: "Scramble table"
    # outputSource: [ScrambleRealigned/table, Scramble/table]
  
  pesr_disc:
    type: File?
    secondaryFiles: [.tbi]
    doc: "PESR discordant reads"
    # outputSource: CollectSVEvidence/disc_out
  
  pesr_disc_index:
    type: File?
    doc: "PESR discordant reads index"
    # outputSource: CollectSVEvidence/disc_out_index
  
  pesr_split:
    type: File?
    secondaryFiles: [.tbi]
    doc: "PESR split reads"
    # outputSource: CollectSVEvidence/split_out
  
  pesr_split_index:
    type: File?
    doc: "PESR split reads index"
    # outputSource: CollectSVEvidence/split_out_index
  
  pesr_sd:
    type: File?
    secondaryFiles: [.tbi]
    doc: "PESR SD information"
    # outputSource: CollectSVEvidence/sd_out
  
  pesr_sd_index:
    type: File?
    doc: "PESR SD information index"
    # outputSource: CollectSVEvidence/sd_out_index
  
  wham_vcf:
    type: File?
    secondaryFiles: [.tbi]
    doc: "Wham VCF"
    # outputSource: Whamg/vcf
  
  wham_index:
    type: File?
    doc: "Wham VCF index"
    # outputSource: Whamg/index
  
  sample_metrics_files:
    type: File[]?
    doc: "Sample metrics files"
    # outputSource: GatherSampleEvidenceMetrics/sample_metrics_files

steps:
  # TODO: Add workflow steps for:
  # 1. LocalizeReads
  # 2. CollectCounts
  # 3. Manta
  # 4. CollectSVEvidence
  # 5. MELT
  # 6. CheckAligner and Scramble
  # 7. RealignSoftClippedReads and ScrambleRealigned if needed
  # 8. Whamg
  # 9. GatherSampleEvidenceMetrics