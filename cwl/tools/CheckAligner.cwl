#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool
label: CheckAligner
doc: "Check if BAM/CRAM was aligned with Dragen 3.7.8"

requirements:
  DockerRequirement:
    dockerPull: $(inputs.gatk_docker)
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}

inputs:
  reads_path:
    type: File
    doc: "BAM or CRAM file"
    secondaryFiles: [.bai, .crai]
    inputBinding:
      prefix: "-I"
  
  reads_index:
    type: File
    doc: "Index file for BAM or CRAM"
    inputBinding:
      prefix: "--read-index"
  
  reference_fasta:
    type: File?
    doc: "Reference FASTA file"
    secondaryFiles: [.fai, ^.dict]
    inputBinding:
      prefix: "-R"
      valueFrom: $(self ? self.path : null)
  
  reference_index:
    type: File?
    doc: "Reference FASTA index file"
  
  reference_dict:
    type: File?
    doc: "Reference FASTA dictionary file"
  
  sample_id:
    type: string
    doc: "Sample ID"
  
  gatk_docker:
    type: string
    doc: "GATK Docker image"

baseCommand: []

arguments:
  - position: 0
    shellQuote: false
    valueFrom: |
      set -euo pipefail

      gatk PrintReadsHeader \
        -I $(inputs.reads_path.path) \
        --read-index $(inputs.reads_index.path) \
        -O $(inputs.sample_id).header.sam \
        $(inputs.reference_fasta ? "-R " + inputs.reference_fasta.path : "")

      awk '$0~"@PG" && $0~"ID: DRAGEN SW build" && $0~"VN: 05.021.604.3.7.8"' $(inputs.sample_id).header.sam \
        | wc -l \
        > is_dragen_3_7_8.txt

outputs:
  header:
    type: File
    outputBinding:
      glob: "$(inputs.sample_id).header.sam"
  
  is_dragen_3_7_8:
    type: int
    outputBinding:
      loadContents: true
      glob: "is_dragen_3_7_8.txt"
      outputEval: $(parseInt(self[0].contents))

hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 1024
    outdirMin: 10240