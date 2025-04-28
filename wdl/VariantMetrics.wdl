version 1.0

import "Structs.wdl"

workflow VariantMetrics {
    input {
        Array[File]? manta_vcfs
        Array[File]? melt_vcfs
        Array[File]? wham_vcfs
        Array[File]? scramble_vcfs
        Array[String] samples
        String prefix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_count
        RuntimeAttr? runtime_attr_merge
    }

    # Process each sample's VCFs to get counts
    scatter (i in range(length(samples))) {
        call CountVariants {
            input:
                sample_id = samples[i],
                manta_vcf = if defined(manta_vcfs) then manta_vcfs[i] else None,
                melt_vcf = if defined(melt_vcfs) then melt_vcfs[i] else None,
                wham_vcf = if defined(wham_vcfs) then wham_vcfs[i] else None,
                scramble_vcf = if defined(scramble_vcfs) then scramble_vcfs[i] else None,
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_count
        }
    }

    # Merge all sample counts into one TSV
    call MergeCounts {
        input:
            count_files = CountVariants.count_file,
            prefix = prefix,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_merge
    }

    output {
        File stratified_counts = MergeCounts.merged_counts
    }
}

task CountVariants {
    input {
        String sample_id
        File? manta_vcf
        File? melt_vcf
        File? wham_vcf
        File? scramble_vcf
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File count_file = "${sample_id}.variant_counts.txt"
    }

    command <<<
        set -euo pipefail

        # Function to process each VCF and output counts
        process_vcf() {
            local vcf=$1
            local caller=$2
            local output=$3

            # Count variants by type and chromosome
            python <<CODE
import gzip

def get_info_field(info_str, field):
    for item in info_str.split(';'):
        if '=' in item and item.split('=')[0] == field:
            return item.split('=')[1]
    return None

counts = {}
with gzip.open('${vcf}', 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        chrom = fields[0]
        info = fields[7]
        svtype = get_info_field(info, 'SVTYPE')
        if svtype:
            key = f"{caller}_{svtype}_{chrom}"
            counts[key] = counts.get(key, 0) + 1

with open(output, 'w') as out:
    for key, count in sorted(counts.items()):
        out.write(f"{key}\t{count}\n")
CODE
        }

        # Initialize output file with header
        echo -e "sample_id\tmetric\tcount" > ~{sample_id}.variant_counts.txt

        # Process each caller's VCF if provided
        ~{if defined(manta_vcf) then "process_vcf ~{manta_vcf} 'manta' 'manta_counts.txt' && cat manta_counts.txt | while read -r metric count; do echo -e '~{sample_id}\t'$metric'\t'$count; done >> ~{sample_id}.variant_counts.txt" else ""}
        ~{if defined(melt_vcf) then "process_vcf ~{melt_vcf} 'melt' 'melt_counts.txt' && cat melt_counts.txt | while read -r metric count; do echo -e '~{sample_id}\t'$metric'\t'$count; done >> ~{sample_id}.variant_counts.txt" else ""}
        ~{if defined(wham_vcf) then "process_vcf ~{wham_vcf} 'wham' 'wham_counts.txt' && cat wham_counts.txt | while read -r metric count; do echo -e '~{sample_id}\t'$metric'\t'$count; done >> ~{sample_id}.variant_counts.txt" else ""}
        ~{if defined(scramble_vcf) then "process_vcf ~{scramble_vcf} 'scramble' 'scramble_counts.txt' && cat scramble_counts.txt | while read -r metric count; do echo -e '~{sample_id}\t'$metric'\t'$count; done >> ~{sample_id}.variant_counts.txt" else ""}
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task MergeCounts {
    input {
        Array[File] count_files
        String prefix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File merged_counts = "${prefix}.stratified_counts.tsv"
    }

    command <<<
        set -euo pipefail

        # Merge all count files and pivot to create final TSV
        python <<CODE
import pandas as pd

# Read all count files
dfs = []
for f in '~{sep="' '" count_files}'.split():
    df = pd.read_csv(f, sep='\t')
    dfs.append(df)

# Combine all dataframes
combined = pd.concat(dfs, ignore_index=True)

# Pivot the table to get metrics as columns
pivoted = combined.pivot(index='sample_id', columns='metric', values='count')
pivoted = pivoted.fillna(0)  # Replace NaN with 0
pivoted = pivoted.astype(int)  # Convert to integers

# Save to TSV
pivoted.to_csv('~{prefix}.stratified_counts.tsv', sep='\t')
CODE
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
} 
