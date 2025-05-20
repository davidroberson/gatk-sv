#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool
label: LocalizeReads
doc: "Move or copy reads from a possibly remote source to a local working directory"

requirements:
  DockerRequirement:
    dockerPull: ubuntu:18.04
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing: []
  ShellCommandRequirement: {}

inputs:
  reads_path:
    type: File
    doc: "BAM or CRAM file"
    inputBinding:
      position: 1
      shellQuote: true
  
  reads_index:
    type: File
    doc: "Index file for BAM or CRAM"
    inputBinding:
      position: 2
      shellQuote: true
  
  move_files:
    type: boolean
    default: false
    doc: "If true, move files instead of copying them"
    inputBinding:
      position: 3
      shellQuote: true

baseCommand: []

arguments:
  - position: 0
    shellQuote: false
    valueFrom: |
      set -exuo pipefail

      if $(inputs.move_files); then
        mv $(inputs.reads_path.path) $(basename "$(inputs.reads_path.basename)")
        mv $(inputs.reads_index.path) $(basename "$(inputs.reads_index.basename)")
      else
        cp $(inputs.reads_path.path) $(basename "$(inputs.reads_path.basename)")
        cp $(inputs.reads_index.path) $(basename "$(inputs.reads_index.basename)")
      fi

outputs:
  output_file:
    type: File
    outputBinding:
      glob: $(basename("$(inputs.reads_path.basename)"))
  
  output_index:
    type: File
    outputBinding:
      glob: $(basename("$(inputs.reads_index.basename)"))

hints:
  ResourceRequirement:
    coresMin: 2
    ramMin: $(3.75 * 1024)
    outdirMin: $(50 + (inputs.move_files ? (inputs.reads_path.size/(1024*1024*1024)) : (inputs.reads_path.size/(1024*1024*1024) * 2)))