# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

GATK-SV is a structural variation discovery pipeline for Illumina short-read whole-genome sequencing (WGS) data. The pipeline is implemented in Workflow Description Language (WDL) and is designed to run on cloud platforms through Cromwell or Terra.

The pipeline has two primary calling modes:
- **Single-sample mode**: Calls SVs on individual samples against a fixed reference panel
- **Joint calling mode**: Performs SV calling on cohorts of samples (recommended for 100+ samples)

The pipeline maximizes SV discovery sensitivity by harmonizing output from multiple callers: Manta, Wham, Scramble, cn.MOPS, and GATK-gCNV.

## Repository Structure

- `/dockerfiles`: Resources for building pipeline docker images
- `/inputs`: Files for generating workflow inputs
  - `/templates`: Input json file templates
  - `/values`: Input values used to populate templates
- `/scripts`: Scripts for testing, building dockers, etc.
- `/src`: Main pipeline scripts organized by module
- `/wdl`: WDL workflow definitions
- `/website`: Documentation website code

## Common Commands

### Python Linting

Run flake8 linting using tox:

```bash
tox -e lint
```

### WDL Validation

Validate WDL files using WOMtool:

```bash
# Install WOMtool
wget -O womtool.jar https://github.com/broadinstitute/cromwell/releases/download/84/womtool-84.jar

# Basic validation of test inputs
./scripts/test/validate.sh -d $(pwd) -j womtool.jar

# Include Terra cohort mode input JSONs validation
./scripts/test/validate.sh -d $(pwd) -j womtool.jar -t
```

Validate WDL files with miniwdl:

```bash
# Install miniwdl
pip install miniwdl

# Validate WDLs
python scripts/test/miniwdl_validation.py --imports-dir wdl wdl/*.wdl
```

### Generate Default Inputs

Generate default inputs for testing:

```bash
scripts/inputs/build_default_inputs.sh -d $(pwd)
```

### Clean Default Inputs

Remove generated default inputs:

```bash
scripts/inputs/clean_default_inputs.sh
```

### Building Docker Images

Build and push Docker images:

```bash
# Install dependencies
pip install termcolor pprint

# Build a specific Docker image
cd scripts/docker
python build_docker.py --targets sv-pipeline --image-tag my-tag --docker-repo myrepo

# Build all Docker images (can take several hours)
python build_docker.py --targets all --image-tag my-tag --docker-repo myrepo
```

Note: Due to license restrictions, you must obtain a license and download `MELTv2.0.5_patch.tar.gz` to `dockerfiles/melt` to build the MELT docker.

### Website Development

To build and run the documentation website locally:

```bash
# Prerequisites
# Install Node.js and Git LFS

# Fetch LFS files
git lfs install
git lfs fetch --all
git lfs pull

# Install dependencies
cd website/
npm install

# Start development server
npm run start

# Build static site
npm run build
```

## Pipeline Architecture

The GATK-SV pipeline consists of multiple modular WDL workflows designed to be run in sequence for joint calling, while single-sample mode is implemented as a single workflow.

### Key Modules (Joint Calling Mode)

1. **Evidence Collection**
   - GatherSampleEvidence: Runs SV discovery tools on each sample
   - GatherBatchEvidence: Merges evidence across samples

2. **Evidence QC** (EvidenceQC)
   - Quality assessment of raw SV evidence

3. **Clustering** (ClusterBatch)
   - Group similar variants across samples

4. **Filtering**
   - FilterBatchSites: Remove false positive variant sites
   - FilterBatchSamples: Identify and filter problematic samples

5. **Genotyping**
   - GenotypeBatch: Re-genotype SVs across all samples
   - RegenotypeCNVs: Specialized genotyping for copy number variants

6. **Complex Variant Resolution**
   - ResolveComplexVariants: Identify and resolve complex SVs
   - RefineComplexVariants: Fine-tune complex variant calls

7. **Annotation** (AnnotateVcf)
   - Add functional annotations to variants

8. **Quality Control** (MainVcfQc)
   - Assess quality of final call set

### Single Sample Mode

Single sample mode uses a simplified workflow that follows similar steps but optimized for individual samples, typically used when comparing against a reference panel.

## Input Requirements

- Docker images for all tools
- Reference files and resources (specified in JSON files in `/inputs/values/`)
- BAM files from aligned sequencing data
- Sample information

## Development Best Practices

When contributing to GATK-SV:

1. Always validate WDLs after making changes using the validation commands above
2. Update Docker images when modifying tool dependencies
3. Follow Python style guidelines (use the lint command to check)
4. When modifying WDLs, test with the provided test inputs
5. Update documentation in the website when adding new features

## Troubleshooting

For common issues and troubleshooting help, refer to the documentation website at https://broadinstitute.github.io/gatk-sv/