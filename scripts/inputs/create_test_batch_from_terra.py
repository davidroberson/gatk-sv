#!/bin/python

# Synopsis:
#  Creates input values for a new test batch from a Terra workspace data table
#

import argparse
import pandas as pd
import json
import sys
import argparse
import io
import json
import os
import sys
import zipfile

import pandas as pd

import firecloud.api as fapi

SAMPLE_KEYS_MAP = {
    "entity:sample_id": "samples",
    "bam_or_cram_file": "bam_or_cram_files",
    "coverage_counts": "counts",
    "manta_vcf": "manta_vcfs",
    "manta_index" : "manta_vcfs_index",
    "pesr_disc": "PE_files",
    "pesr_disc_index": "PE_files_index",
    "pesr_sd": "SD_files",
    "pesr_sd_index": "SD_files_index",
    "pesr_split": "SR_files",
    "pesr_split_index": "SR_files_index",
    "scramble_vcf": "scramble_vcfs",
    "scramble_index": "scramble_vcfs_index",
    "wham_vcf": "wham_vcfs",
    "wham_index": "wham_vcfs_index"
}

SAMPLE_SET_KEYS_MAP = {
    "entity:sample_set_id": "name",
    "clustered_depth_vcf": "merged_depth_vcf",
    "clustered_depth_vcf_index": "merged_depth_vcf_index",
    "clustered_manta_vcf": "merged_manta_vcf",
    "clustered_manta_vcf_index": "merged_manta_vcf_index",
    "clustered_scramble_vcf": "merged_scramble_vcf",
    "clustered_scramble_vcf_index": "merged_scramble_vcf_index",
    "clustered_wham_vcf": "merged_wham_vcf",
    "clustered_wham_vcf_index": "merged_wham_vcf_index",
    "contig_ploidy_model_tar": "contig_ploidy_model_tar",
    "cutoffs": "cutoffs",
    "genotyped_depth_vcf": "genotyped_depth_vcf",
    "genotyped_pesr_vcf": "genotyped_pesr_vcf",
    "median_cov": "medianfile",
    "merged_BAF": "merged_baf_file",
    "merged_BAF_index": "merged_baf_file_index",
    "merged_PE": "merged_disc_file",
    "merged_PE_index": "merged_disc_file_index",
    "merged_SR": "merged_split_file",
    "merged_SR_index": "merged_split_file_index",
    "merged_bincov": "merged_coverage_file",
    "merged_bincov_index": "merged_coverage_file_index",
    "merged_dels": "del_bed",
    "merged_dups": "dup_bed",
    "metrics": "evidence_metrics",
    "metrics_common": "evidence_metrics_common",
    "outlier_filtered_depth_vcf": "filtered_depth_vcf",
    "outlier_filtered_depth_vcf_index": "filtered_depth_vcf_index",
    "outlier_filtered_pesr_vcf": "filtered_pesr_vcf",
    "outlier_filtered_pesr_vcf_index": "filtered_pesr_vcf_index",
    "regeno_coverage_medians": "regeno_coverage_medians",
    "sites_filtered_depth_vcf": "sites_filtered_depth_vcf",
    "sites_filtered_manta_vcf": "sites_filtered_manta_vcf",
    "sites_filtered_scramble_vcf": "sites_filtered_scramble_vcf",
    "sites_filtered_wham_vcf": "sites_filtered_wham_vcf",
    "sr_background_fail": "cluster_background_fail_lists",
    "sr_bothside_pass": "cluster_bothside_pass_lists",
    "std_manta_vcf_tar": "std_manta_vcf_tar",
    "std_scramble_vcf_tar": "std_scramble_vcf_tar",
    "std_wham_vcf_tar": "std_wham_vcf_tar",
    "trained_PE_metrics": "PE_metrics",
    "trained_SR_metrics": "SR_metrics",
    "trained_genotype_depth_depth_sepcutoff": "genotype_depth_depth_sepcutoff",
    "trained_genotype_depth_pesr_sepcutoff": "genotype_depth_pesr_sepcutoff",
    "trained_genotype_pesr_depth_sepcutoff": "genotype_pesr_depth_sepcutoff",
    "trained_genotype_pesr_pesr_sepcutoff": "genotype_pesr_pesr_sepcutoff"
}

WORKSPACE_DATA_KEY_MAP = {
    "cohort_ped_file": "ped_file"
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples-tsv", help="Sample data table tsv")
    parser.add_argument("--sample-set-entity-tsv", help="Sample set data table entity tsv")
    parser.add_argument("--workspace-tsv", help="Workspace data tsv")
    parser.add_argument("--new-bucket-name", help="If used, changes bucket name for all files")
    args = parser.parse_args()

    df_samples = pd.read_csv(args.samples_tsv, sep='\t')\
        .rename(columns={"entity:sample_id": "sample_id"})\
        .set_index('sample_id')
    print(df_samples)
    df_entity = pd.read_csv(args.sample_set_entity_tsv, sep='\t')\
        .rename(columns={"entity:sample_set_id": "sample_set_id"})\
        .set_index('sample_set_id')
    print(df_entity)
    if df_entity.shape[0] != 1:
        raise ValueError(f"Expected exactly 1 sample set in entity tsv but found {df_entity.shape[0]}")
    df_workspace = pd.read_csv(args.workspace_tsv, sep='\t')
    df_workspace = df_workspace.rename(columns={column: column.replace("workspace:", "") for column in
                                                df_workspace.columns if column.startswith("workspace:")})
    print(df_workspace)


if __name__ == "__main__":
    main()
