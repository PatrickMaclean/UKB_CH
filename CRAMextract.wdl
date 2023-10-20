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
        File hla_out = cram_munging.mhc_cram
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
samtools view -b -T reference.fasta -L regions.bed -o output.bam file.cram && samtools index output.bam

samtools view --reference GRCh38_full_analysis_set_plus_decoy_hla.fa -L CHIP_genes_38.bed  -o 1090728_23143_0_0_CHIPregions.bam 1090728_23143_0_0.cram
samtools index 1090728_23143_0_0_CHIPregions.bam

        bcftools filter --regions chr1:159612289-159814589 ~{vcf_file} -Ou | \
            bcftools view -f 'FILTER=PASS' -Ou | \
            bcftools norm -m -any --check-ref w -f ~{ref_genome} -Ou | \
            bcftools plugin fill-tags -Ou | \
            bcftools view --max-af 0.01:minor -Ou | \
            bcftools annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -Oz > ukb_wgs_mhc.vcf.gz
    >>>

    runtime {
        docker: "monsieurbl/bcftools"
        dx_instance_type: "mem1_ssd1_v2_x8"
    }

    output {
        File exome_bam = "
        File exome_bai = 
        File mhc_vcf = "ukb_wgs_mhc.vcf.gz"
    }

}
