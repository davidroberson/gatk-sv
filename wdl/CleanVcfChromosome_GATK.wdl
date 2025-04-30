version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "FormatVcfForGatk.wdl" as fvcf

workflow CleanVcfChromosome_GATK {
  input {
    File vcf
    String contig
    File ped_file
    String prefix

    # Required parameters for ploidy handling
    String chr_x
    String chr_y
    File? male_samples_list

    # Reference files for mobile element deletion rescue
    File? LINE1_reference
    File? HERVK_reference

    # Runtime parameters
    String sv_pipeline_docker
    String sv_base_mini_docker
    String gatk_docker

    # Runtime attributes for each task
    RuntimeAttr? runtime_attr_create_ploidy_table
    RuntimeAttr? runtime_attr_format_vcf_for_gatk
    RuntimeAttr? runtime_attr_preprocess_vcf
    RuntimeAttr? runtime_attr_overlapping_cnvs
    RuntimeAttr? runtime_attr_multiallelic_cnvs
    RuntimeAttr? runtime_attr_abnormal_allosomes
    RuntimeAttr? runtime_attr_overlapping_multiallelics
    RuntimeAttr? runtime_attr_postprocess_vcf
    
    # Runtime attributes for additional tasks
    RuntimeAttr? runtime_attr_drop_redundant_cnvs
    RuntimeAttr? runtime_override_sort_drop_redundant_cnvs
    RuntimeAttr? runtime_attr_stitch_fragmented_cnvs
    RuntimeAttr? runtime_attr_rescue_me_dels
    RuntimeAttr? runtime_attr_add_high_fp_rate_filters
    RuntimeAttr? runtime_attr_final_cleanup
    RuntimeAttr? runtime_attr_format
    
    # Other parameters
    Boolean use_hail = false
    String? gcs_project
  }

  # Create ploidy table from PED file
  call CreatePloidyTable {
    input:
      ped_file = ped_file,
      chr_x = chr_x,
      chr_y = chr_y,
      prefix = "~{prefix}.ploidy",
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_create_ploidy_table
  }

  # Format VCF for GATK
  call FormatVcfForGatk {
    input:
      vcf = vcf,
      prefix = "~{prefix}.format_for_gatk",
      ploidy_table = CreatePloidyTable.out,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_format_vcf_for_gatk
  }

  # Preprocess VCF
  call PreprocessVcf {
    input:
      vcf = FormatVcfForGatk.out_vcf,
      prefix = "~{prefix}.preprocess",
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_preprocess_vcf
  }

  # Run the GATK Overlapping CNVs revision
  call GatkSVReviseOverlappingCnvs {
    input:
      vcf = PreprocessVcf.out_vcf,
      prefix = "~{prefix}.overlapping_cnvs",
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_overlapping_cnvs
  }

  # Run the GATK Multiallelic CNVs revision
  call GatkSVReviseMultiallelicCnvs {
    input:
      vcf = GatkSVReviseOverlappingCnvs.out_vcf,
      prefix = "~{prefix}.multiallelic_cnvs",
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_multiallelic_cnvs
  }

  # Run the GATK Abnormal Allosomes revision
  call GatkSVReviseAbnormalAllosomes {
    input:
      vcf = GatkSVReviseMultiallelicCnvs.out_vcf,
      prefix = "~{prefix}.abnormal_allosomes",
      male_samples_list = male_samples_list,
      chr_x = chr_x,
      chr_y = chr_y,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_abnormal_allosomes
  }

  # Run the GATK Overlapping Multiallelics revision
  call GatkSVReviseOverlappingMultiallelics {
    input:
      vcf = GatkSVReviseAbnormalAllosomes.out_vcf,
      prefix = "~{prefix}.overlapping_multiallelics",
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_overlapping_multiallelics
  }

  # Postprocess the VCF
  call PostprocessVcf {
    input:
      vcf = GatkSVReviseOverlappingMultiallelics.out_vcf,
      prefix = "~{prefix}.postprocess",
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_postprocess_vcf
  }

  # Add the remaining tasks from the original workflow

  # Drop redundant CNVs
  call DropRedundantCnvs {
    input:
      vcf = PostprocessVcf.out_vcf,
      prefix = "~{prefix}.drop_redundant_cnvs",
      contig = contig,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_drop_redundant_cnvs
  }

  # Sort the VCF
  call MiniTasks.SortVcf as SortDropRedundantCnvs {
    input:
      vcf = DropRedundantCnvs.out,
      outfile_prefix = "~{prefix}.drop_redundant_cnvs.sorted",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_override_sort_drop_redundant_cnvs
  }

  # Stitch fragmented CNVs
  call StitchFragmentedCnvs {
    input:
      vcf = SortDropRedundantCnvs.out,
      prefix = "~{prefix}.stitch_fragmented_cnvs",
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_stitch_fragmented_cnvs
  }

  # Conditional call to RescueMobileElementDeletions if references are provided
  if (defined(LINE1_reference) && defined(HERVK_reference)) {
    call RescueMobileElementDeletions {
      input:
        vcf = StitchFragmentedCnvs.stitched_vcf_shard,
        prefix = "~{prefix}.rescue_me_dels",
        LINE1 = select_first([LINE1_reference]),
        HERVK = select_first([HERVK_reference]),
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_rescue_me_dels
    }
  }

  # Use the output from RescueMobileElementDeletions if it ran, otherwise use the output from StitchFragmentedCnvs
  File rescue_me_out = select_first([RescueMobileElementDeletions.out, StitchFragmentedCnvs.stitched_vcf_shard])

  # Add high FDR filters
  call AddHighFDRFilters {
    input:
      vcf = rescue_me_out,
      prefix = "~{prefix}.high_fdr_filtered",
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_add_high_fp_rate_filters
  }

  # Final cleanup
  call FinalCleanup {
    input:
      vcf = AddHighFDRFilters.out,
      contig = contig,
      prefix = "~{prefix}.final_cleanup",
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_final_cleanup
  }

  # Format the final VCF
  call fvcf.FormatVcf {
    input:
      vcf = FinalCleanup.final_cleaned_shard,
      ploidy_table = CreatePloidyTable.out,
      output_prefix = "~{prefix}.final_format",
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_format
  }
  
  output {
    File out = FormatVcf.out
    File out_idx = FormatVcf.out_index
  }
}

# Remove CNVs that are redundant with CPX events or other CNVs
task DropRedundantCnvs {
  input {
    File vcf
    String prefix
    String contig
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(vcf, "GiB")
  Int cpu_cores = 2
  RuntimeAttr default_attr = object {
    mem_gb: 3.75 + input_size * 1.5,
    disk_gb: ceil(100.0 + input_size * 2.0),
    cpu_cores: cpu_cores,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  
  command <<<
    set -euo pipefail
    /opt/sv-pipeline/04_variant_resolution/scripts/resolve_cpx_cnv_redundancies.py \
      ~{vcf} ~{prefix}.vcf.gz --temp-dir ./tmp
  >>>
  
  output {
    File out = "~{prefix}.vcf.gz"
  }
  
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

# Stitch fragmented RD-only calls found in 100% of the same samples
task StitchFragmentedCnvs {
  input {
    File vcf
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(vcf, "GB")
  RuntimeAttr default_attr = object {
    mem_gb: 7.5,
    disk_gb: ceil(10.0 + input_size * 2),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

  command <<<
    set -euo pipefail
    echo "First pass..."
    java -Xmx~{java_mem_mb}M -jar ${STITCH_JAR} 0.2 200000 0.2 ~{vcf} \
      | bgzip \
      > tmp.vcf.gz
    rm ~{vcf}
    echo "Second pass..."
    java -Xmx~{java_mem_mb}M -jar ${STITCH_JAR} 0.2 200000 0.2 tmp.vcf.gz \
      | bgzip \
      > ~{prefix}.vcf.gz
  >>>

  output {
    File stitched_vcf_shard = "~{prefix}.vcf.gz"
  }
  
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

# Rescue mobile element deletions
task RescueMobileElementDeletions {
  input {
    File vcf
    String prefix
    File LINE1
    File HERVK
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(vcf, "GiB")
  RuntimeAttr default_attr = object {
    mem_gb: 3.75 + input_size * 1.5,
    disk_gb: ceil(100.0 + input_size * 3.0),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    python <<CODE
import os
import pysam
fin=pysam.VariantFile("~{vcf}")
fo=pysam.VariantFile("~{prefix}.bnd_del.vcf.gz", 'w', header = fin.header)
for record in fin:
    if record.info['SVTYPE'] in ['BND'] and record.info['STRANDS']=="+-" and record.chrom == record.info['CHR2'] and record.info['END2'] - record.start < 10000:
        record.info['SVLEN'] = record.info['END2'] - record.start
        fo.write(record)
fin.close()
fo.close()
CODE

    tabix -p vcf ~{prefix}.bnd_del.vcf.gz

    svtk vcf2bed ~{prefix}.bnd_del.vcf.gz -i ALL --include-filters ~{prefix}.bnd_del.bed
    bgzip ~{prefix}.bnd_del.bed

    bedtools coverage -wo -a ~{prefix}.bnd_del.bed.gz -b ~{LINE1} | awk '{if ($NF>.5) print}' | cut -f4 | sed -e 's/$/\tDEL\tPASS\toverlap_LINE1/' > manual_revise.MEI_DEL_from_BND.SVID_SVTYPE_FILTER_INFO.tsv
    bedtools coverage -wo -a ~{prefix}.bnd_del.bed.gz -b ~{HERVK} | awk '{if ($NF>.5) print}' | cut -f4 | sed -e 's/$/\tDEL\tPASS\toverlap_HERVK/' >> manual_revise.MEI_DEL_from_BND.SVID_SVTYPE_FILTER_INFO.tsv

    python <<CODE
import pysam
def SVID_MEI_DEL_readin(MEI_DEL_reset):
    out={}
    fin=open(MEI_DEL_reset)
    for line in fin:
        pin=line.strip().split()
        if not pin[0] in out.keys():
            out[pin[0]] = pin[3]
    fin.close()
    return out

hash_MEI_DEL_reset = SVID_MEI_DEL_readin("manual_revise.MEI_DEL_from_BND.SVID_SVTYPE_FILTER_INFO.tsv")
fin=pysam.VariantFile("~{vcf}")
fo=pysam.VariantFile("~{prefix}.vcf.gz", 'w', header = fin.header)
for record in fin:
    if record.id in hash_MEI_DEL_reset.keys():
        del record.filter['UNRESOLVED']
        record.info['SVTYPE'] = 'DEL'
        record.info['SVLEN'] = record.info['END2'] - record.start
        record.stop = record.info['END2']
        record.info.pop("CHR2")
        record.info.pop("END2")
        record.info.pop("UNRESOLVED_TYPE")
        if hash_MEI_DEL_reset[record.id] == 'overlap_LINE1':
            record.alts = ('<DEL:ME:LINE1>',)
        if hash_MEI_DEL_reset[record.id] == 'overlap_HERVK':
            record.alts = ('<DEL:ME:HERVK>',)
    fo.write(record)
fin.close()
fo.close()
CODE
  >>>

  output {
    File out = "~{prefix}.vcf.gz"
  }
  
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

# Add FILTER status for pockets of variants with high FP rate
task AddHighFDRFilters {
  input {
    File vcf
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(vcf, "GiB")
  RuntimeAttr default_attr = object {
    mem_gb: 3.75,
    disk_gb: ceil(10.0 + input_size * 3.0),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    python <<CODE
import pysam
with pysam.VariantFile("~{vcf}", 'r') as fin:
  header = fin.header
  header.add_line("##FILTER=<ID=HIGH_ALGORITHM_FDR,Description=\"Categories of variants with low precision including Wham-only deletions and certain Scramble SVAs\">")
  with pysam.VariantFile("~{prefix}.vcf.gz", 'w', header=header) as fo:
    for record in fin:
        if (record.info['ALGORITHMS'] == ('wham',) and record.info['SVTYPE'] == 'DEL') or \
          (record.info['ALGORITHMS'] == ('scramble',) and record.info['HIGH_SR_BACKGROUND'] and record.alts == ('<INS:ME:SVA>',)):
            record.filter.add('HIGH_ALGORITHM_FDR')
        fo.write(record)
CODE
  >>>

  output {
    File out = "~{prefix}.vcf.gz"
  }
  
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

# Final VCF cleanup
task FinalCleanup {
  input {
    File vcf
    String contig
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(vcf, "GB")
  Float base_disk_gb = 10.0
  Float base_mem_gb = 2.0
  Float input_mem_scale = 3.0
  Float input_disk_scale = 5.0
  RuntimeAttr default_attr = object {
    mem_gb: base_mem_gb + input_size * input_mem_scale,
    disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -eu -o pipefail
    
    /opt/sv-pipeline/04_variant_resolution/scripts/rename_after_vcfcluster.py \
      --chrom ~{contig} \
      --prefix ~{prefix} \
      ~{vcf} stdout \
      | bcftools annotate --no-version -e 'SVTYPE=="CNV" && SVLEN<5000' -x INFO/MEMBERS -Oz -o ~{prefix}.vcf.gz
    tabix ~{prefix}.vcf.gz
  >>>

  output {
    File final_cleaned_shard = "~{prefix}.vcf.gz"
    File final_cleaned_shard_idx = "~{prefix}.vcf.gz.tbi"
  }
  
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

# Create ploidy table from PED file
task CreatePloidyTable {
  input {
    File ped_file
    String chr_x
    String chr_y
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 2.0,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    
    # Create ploidy table from PED file
    python <<CODE
    import pandas as pd
    
    # Read PED file
    ped = pd.read_csv("~{ped_file}", sep='\t', header=None)
    ped.columns = ['fam_id', 'sample_id', 'father_id', 'mother_id', 'sex', 'phenotype']
    
    # Map sex to ploidy for X and Y chromosomes
    # In PED format: 1=male, 2=female, other=unknown
    output_lines = []
    
    for idx, row in ped.iterrows():
        sample = row['sample_id']
        sex = row['sex']
        
        # Set ploidy for X chromosome
        if sex == 1:  # Male
            output_lines.append(f"{sample}\t~{chr_x}\t1")
        elif sex == 2:  # Female
            output_lines.append(f"{sample}\t~{chr_x}\t2")
        else:  # Unknown, assume diploid
            output_lines.append(f"{sample}\t~{chr_x}\t2")
        
        # Set ploidy for Y chromosome
        if sex == 1:  # Male
            output_lines.append(f"{sample}\t~{chr_y}\t1")
        elif sex == 2:  # Female
            output_lines.append(f"{sample}\t~{chr_y}\t0")
        else:  # Unknown, assume zero copies for safety
            output_lines.append(f"{sample}\t~{chr_y}\t0")
        
        # Set ploidy for all other chromosomes as 2
        output_lines.append(f"{sample}\tALL\t2")
    
    # Write to output file
    with open("~{prefix}.txt", "w") as f:
        f.write("\n".join(output_lines))
    CODE
  >>>

  output {
    File out = "~{prefix}.txt"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

# Format VCF for GATK processing
task FormatVcfForGatk {
  input {
    File vcf
    File ploidy_table
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: ceil(size(vcf, "GB") * 3 + 10),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    
    # Format VCF for GATK processing
    python /opt/sv-pipeline/scripts/format_svtk_vcf_for_gatk.py \
      --vcf ~{vcf} \
      --out ~{prefix}.vcf.gz \
      --ploidy-table ~{ploidy_table}
    
    tabix -p vcf ~{prefix}.vcf.gz
  >>>

  output {
    File out_vcf = "~{prefix}.vcf.gz"
    File out_vcf_idx = "~{prefix}.vcf.gz.tbi"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: "~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_attr.disk_gb, default_attr.disk_gb])} HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
