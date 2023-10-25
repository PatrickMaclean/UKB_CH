version 1.0

## notes - https://github.com/openwdl/learn-wdl/blob/master/3_genomic_tool_pipelines/3_haplotype_caller_GATK/HaplotypeCaller.aws.wdl

# WORKFLOW DEFINITION

workflow CRAMextract {
    input {
        File cram_file
        File cram_index_file
        File CHIP_regions       #CHIP_genes_38.bed
        File ref_genome         #GRCh38_full_analysis_set_plus_decoy_hla.fa.gz
        File ref_genome_index   #GRCh38_full_analysis_set_plus_decoy_hla.fa.gz.fai
        File ref_genome_gzi     #GRCh38_full_analysis_set_plus_decoy_hla.fa.gz.gzi
    }

String input_basename = basename(input_cram, ".cram")

    call cram_munging { 
        input: cram_file = cram_file,
        cram_index_file = cram_file+".crai",
        CHIP_regions = CHIP_regions,
        ref_genome = ref_genome,
        ref_genome_index = ref_genome_index,
        ref_genome_gzi = ref_genome_gzi,
        sample_name = input_basename
    }

    output {
        File exome_CHIPregions_bam = cram_munging."~{sample_name}_CHIPregions.bam"
        File exome_CHIPregions_bam_bai = "output.bam.bai"
    }
  
}

# TASK DEFINITIONS

task cram_munging {
    input {
        File cram_file
        File cram_index_file
        File CHIP_regions       #CHIP_genes_38.bed
        File ref_genome         #GRCh38_full_analysis_set_plus_decoy_hla.fa.gz
        File ref_genome_index   #GRCh38_full_analysis_set_plus_decoy_hla.fa.gz.fai
        File ref_genome_gzi     #GRCh38_full_analysis_set_plus_decoy_hla.fa.gz.gzi
    }

    command <<<
samtools view --reference ~{ref_genome} \
  -L ~{CHIP_regions} \
  -o "~{sample_name}_CHIPregions.bam" \
  ~{cram_file} && \
  samtools index "~{sample_name}_CHIPregions.bam"
    >>>

    runtime {
        docker: "monsieurbl/bcftools"
        dx_instance_type: "mem1_ssd1_v2_x8"
    }

    output {
        File exome_bam = "~{sample_name}_CHIPregions.bam"
        File exome_bai = "~{sample_id}_CHIPregions.bam.bai"
    }

}
