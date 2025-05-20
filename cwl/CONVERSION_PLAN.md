# GATK-SV WDL to CWL Conversion Plan

## Overview

This document outlines the plan for converting the GATK-SV pipeline from Workflow Description Language (WDL) to Common Workflow Language (CWL). We'll start with the GatherSampleEvidence workflow, which is a key component in the GATK-SV pipeline.

## Phase 1: GatherSampleEvidence Workflow Conversion

### 1. Directory Structure

```
cwl/
├── GatherSampleEvidence.cwl            # Main workflow
├── Structs.yml                         # Type definitions (equivalent to Structs.wdl)
├── tools/                              # Tasks from main workflow
│   ├── LocalizeReads.cwl
│   ├── CheckAligner.cwl
│   └── RealignSoftClippedReads.cwl
└── workflows/                          # Imported sub-workflows
    ├── CollectCoverage.cwl
    ├── CollectSVEvidence.cwl
    ├── Manta.cwl
    ├── MELT.cwl
    ├── Scramble.cwl
    ├── Whamg.cwl
    └── GatherSampleEvidenceMetrics.cwl
```

### 2. Conversion Steps

1. **Create Type Definitions**:
   - Define runtime attributes
   - Define common input/output types
   - Create Structs.yml for shared data types

2. **Convert Individual Tasks**:
   - LocalizeReads
   - CheckAligner 
   - RealignSoftClippedReads

3. **Convert Import Workflows**:
   - Start with simplest workflows (CollectCoverage, Manta)
   - Move to more complex workflows (MELT, Scramble)
   - Complete remaining workflows

4. **Create Main Workflow**:
   - Define workflow inputs/outputs
   - Add steps for each task/sub-workflow
   - Implement conditional execution logic
   - Connect inputs and outputs between steps

5. **Testing**:
   - Create test data set
   - Validate syntax with cwltool
   - Test execution with small test dataset
   - Compare outputs with WDL version

### 3. CWL-Specific Implementation Notes

#### Handling WDL-Specific Features

1. **Optional Inputs and Outputs**:
   - Use `null` for optional values in CWL
   - Include appropriate type definitions with `?` for optional types

2. **Conditional Execution**:
   - Use `when` expressions in CWL steps
   - Implement scatter-gather patterns where needed

3. **Select First Logic**:
   - Replace WDL's `select_first()` with default values in CWL inputs

4. **Runtime Attributes**:
   - Define resource requirements directly in CWL
   - Create reusable requirement sets for common patterns

### 4. Detailed Implementation Plan

#### A. GatherSampleEvidence.cwl (Main Workflow)

1. Define workflow-level inputs/outputs
2. Implement conditional logic for running different SV callers
3. Connect sub-workflows and tasks
4. Handle file type determination (BAM vs CRAM)

#### B. Tools (Tasks)

1. **LocalizeReads.cwl**:
   - Handle file localization
   - Implement file movement vs copying option

2. **CheckAligner.cwl**:
   - Parse BAM headers for aligner information
   - Detect Dragen 3.7.8

3. **RealignSoftClippedReads.cwl**:
   - Extract and realign soft-clipped reads
   - Coordinate with BWA for realignment

#### C. Workflows (Sub-workflows)

1. **CollectCoverage.cwl**:
   - Process read counts over intervals
   - Handle memory requirements

2. **Manta.cwl**, **MELT.cwl**, **Scramble.cwl**, **Whamg.cwl**:
   - Implement SV callers
   - Handle caller-specific parameters and outputs

3. **CollectSVEvidence.cwl**:
   - Collect paired-end and split-read evidence
   - Process discordant read pairs

4. **GatherSampleEvidenceMetrics.cwl**:
   - Generate QC metrics
   - Format output for consumption by later steps

## Phase 2: Single Sample Pipeline Conversion

After successfully converting and testing the GatherSampleEvidence workflow, we'll continue with the remaining components of the single sample mode pipeline in this order:

1. EvidenceQC
2. GatherBatchEvidence
3. ClusterBatch 
4. FilterBatch
5. GenotypeBatch
6. MakeCohortVcf
7. AnnotateVcf
8. Main GATKSVPipelineSingleSample workflow

## Phase 3: Full Pipeline Conversion

After the single sample mode is working, we'll convert the remaining workflows for joint calling mode.

## Testing Strategy

1. **Unit Testing**:
   - Test each task and sub-workflow individually
   - Verify outputs match WDL version

2. **Integration Testing**:
   - Test complete GatherSampleEvidence workflow
   - Compare results against WDL outputs

3. **Validation**:
   - Use cwltool --validate to ensure CWL is well-formed
   - Run with --print-rdf to inspect workflow structure

## Success Criteria

1. CWL workflow passes validation
2. Execution completes successfully
3. Outputs match WDL workflow outputs within tolerance
4. Documentation is complete
5. Test suite is implemented