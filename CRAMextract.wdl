version 1.0

workflow CRAMextract {
    input {
        File cram_file
        File cram_index_file
        File CHIP_regions       #CHIP_genes_38.bed
        File ref_genome         #GRCh38_full_analysis_set_plus_decoy_hla.fa.gz
        File ref_genome_index   #GRCh38_full_analysis_set_plus_decoy_hla.fa.gz.fai
        File ref_genome_gzi     #GRCh38_full_analysis_set_plus_decoy_hla.fa.gz.gzi
    }

    call cram_munging { 
        input: cram_file = cram_file, index_file = index_file, ref_genome = ref_genome, ref_genome_index = ref_genome_index, ref_genome_gzi = ref_genome_gzi
    }

    output {
        File bam_out = cram_munging.exome_bam
    }
  
}


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
  -o output.bam \
  ~{cram_file} && \
  samtools index output.bam
    >>>

    runtime {
        docker: "monsieurbl/bcftools"
        dx_instance_type: "mem1_ssd1_v2_x8"
    }

    output {
        File exome_bam = "output.bam"
        File exome_bai = "output.bam.bai"
    }

}
