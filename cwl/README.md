# GATK-SV CWL Pipeline

This directory contains a Common Workflow Language (CWL) implementation of the GATK-SV pipeline.

## Status

This is an experimental conversion from WDL to CWL, starting with the GatherSampleEvidence workflow. This work is in progress and not yet ready for production use.

## Directory Structure

```
cwl/
├── GatherSampleEvidence.cwl        # Main workflow (in progress)
├── Structs.yml                     # Type definitions
├── CONVERSION_PLAN.md              # Detailed conversion plan
├── tools/                          # Individual tasks
│   ├── LocalizeReads.cwl
│   ├── CheckAligner.cwl
│   └── RealignSoftClippedReads.cwl
└── workflows/                      # Sub-workflows
    ├── CollectCoverage.cwl         # (Coming soon)
    ├── CollectSVEvidence.cwl       # (Coming soon)
    ├── Manta.cwl                   # (Coming soon)
    ├── MELT.cwl                    # (Coming soon)
    ├── Scramble.cwl                # (Coming soon)
    ├── Whamg.cwl                   # (Coming soon)
    └── GatherSampleEvidenceMetrics.cwl # (Coming soon)
```

## Testing

This CWL conversion is still in the early stages and does not yet have testing infrastructure.

## Usage

Once completed, the CWL pipeline will offer a cross-platform alternative to the WDL implementation, allowing users to run the pipeline on platforms that support CWL.

## Next Steps

1. Complete the GatherSampleEvidence.cwl workflow
2. Implement and test all sub-workflows
3. Create test inputs and validation process
4. Extend to other components of the GATK-SV pipeline

## References

- [CWL Specification](https://www.commonwl.org/v1.2/)
- [Original GATK-SV WDL Pipeline](https://github.com/broadinstitute/gatk-sv)