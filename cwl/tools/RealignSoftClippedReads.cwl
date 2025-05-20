#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool
label: RealignSoftClippedReads
doc: "Extract and realign soft-clipped reads from BAM/CRAM file"

requirements:
  DockerRequirement:
    dockerPull: $(inputs.sv_base_mini_docker)
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.reads_path)
      - $(inputs.reads_index)
      - $(inputs.scramble_table)
      - $(inputs.reference_fasta)
      - $(inputs.reference_index)
      - $(inputs.reference_bwa_alt)
      - $(inputs.reference_bwa_amb)
      - $(inputs.reference_bwa_ann)
      - $(inputs.reference_bwa_bwt)
      - $(inputs.reference_bwa_pac)
      - $(inputs.reference_bwa_sa)

inputs:
  reads_path:
    type: File
    doc: "BAM or CRAM file"
    secondaryFiles: [.bai, .crai]
  
  reads_index:
    type: File
    doc: "Index file for BAM or CRAM"
  
  scramble_table:
    type: File
    doc: "Scramble table file from previous Scramble step"
  
  is_bam:
    type: boolean
    doc: "Whether input is BAM (true) or CRAM (false)"
  
  sample_id:
    type: string
    doc: "Sample ID"
  
  reference_fasta:
    type: File
    doc: "Reference FASTA file"
    secondaryFiles: [.fai]
  
  reference_index:
    type: File
    doc: "Reference FASTA index file"
  
  reference_bwa_alt:
    type: File
    doc: "BWA ALT file"
  
  reference_bwa_amb:
    type: File
    doc: "BWA AMB file"
  
  reference_bwa_ann:
    type: File
    doc: "BWA ANN file"
  
  reference_bwa_bwt:
    type: File
    doc: "BWA BWT file"
  
  reference_bwa_pac:
    type: File
    doc: "BWA PAC file"
  
  reference_bwa_sa:
    type: File
    doc: "BWA SA file"
  
  sv_base_mini_docker:
    type: string
    doc: "SV base mini Docker image"
  
  use_ssd:
    type: boolean
    default: true
    doc: "Whether to use SSD storage"

baseCommand: []

arguments:
  - position: 0
    shellQuote: false
    valueFrom: |
      set -exuo pipefail
      
      # Get insertion intervals
      zcat $(inputs.scramble_table.path) \
        | sed 1d \
        | cut -f1 \
        | tr ':' '\t' \
        | awk -F'\t' -v OFS='\t' '{print $1,$2,$2+1}' \
        | sort -k1,1V -k2,2n \
        | bedtools slop -i - -g $(inputs.reference_index.path) -b 150 \
        | bedtools merge \
        > intervals.bed
        
      mkdir -p tmp/
      samtools view --header-only $(inputs.reads_path.path) > header.sam
      
      N_CORES=$(nproc)
      
      time samtools view --no-header \
        -T $(inputs.reference_fasta.path) \
        -ML intervals.bed \
        $(inputs.reads_path.path) \
        | awk -F'\t' -v OFS='\t' '$6~"S"' \
        | sort -u \
        | cat header.sam - \
        | samtools fastq \
        > reads.fastq
        
      bwa mem -H header.sam -K 100000000 -v 3 -t ${N_CORES} -Y $(inputs.reference_fasta.path) reads.fastq \
        | samtools sort -T tmp \
        | samtools view -1 -h -O BAM -o $(inputs.sample_id).realign_soft_clipped_reads.bam
        
      samtools index -@${N_CORES} $(inputs.sample_id).realign_soft_clipped_reads.bam

outputs:
  out:
    type: File
    secondaryFiles: [.bai]
    outputBinding:
      glob: "$(inputs.sample_id).realign_soft_clipped_reads.bam"
  
  out_index:
    type: File
    outputBinding:
      glob: "$(inputs.sample_id).realign_soft_clipped_reads.bam.bai"

hints:
  ResourceRequirement:
    coresMin: 4
    ramMin: $(20 * 1024)
    outdirMin: $(10240 + inputs.reads_path.size * 2 + inputs.reference_bwa_bwt.size + inputs.reference_bwa_pac.size + inputs.reference_bwa_sa.size)
    tmpdirMin: 10240