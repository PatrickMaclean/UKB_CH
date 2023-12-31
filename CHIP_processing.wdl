version 1.0

#need to change docker image in the hla_calling task runtime section

workflow CHIP_processing {
    input {
        File cram_file
        File cram_file_index 
        File CHIP_regions
        File ref_genome         #https://github.com/broadinstitute/gatk/raw/master/src/test/resources/large/Homo_sapiens_assembly38.fasta.gz
        File ref_genome_index   #https://github.com/broadinstitute/gatk/raw/master/src/test/resources/large/Homo_sapiens_assembly38.fasta.gz.fai
        File ref_genome_dict    #https://github.com/broadinstitute/gatk/raw/master/src/test/resources/large/Homo_sapiens_assembly38.dict
        File ref_genome_gzi     #https://github.com/broadinstitute/gatk/raw/master/src/test/resources/large/Homo_sapiens_assembly38.fasta.gz.gzi
    }

    call cram_extract_index { 
        input: cram_file = cram_file, cram_file_index = cram_file_index, CHIP_regions = CHIP_regions, ref_genome = ref_genome, ref_genome_index = ref_genome_index, ref_genome_dict = ref_genome_dict, ref_genome_gzi = ref_genome_gzi
    }


    output {
        File extracted_regions = cram_extract_index.extracted_regions
        File extracted_regions_index = cram_extract_index.extracted_regions_index
    }
  
}

task cram_extract_index {
    
    input {
        File cram_file
        File cram_file_index
        File CHIP_regions
        File ref_genome
        File ref_genome_index
        File ref_genome_dict
        File ref_genome_gzi
    }

    command <<<
        set -x -e -o pipefail
        time samtools view --reference ~{ref_genome} -L ~{CHIP_regions} -o extracted_regions.bam ~{cram_file} 
        time samtools index extracted_regions.bam
    >>>

    runtime {
        docker: "dx://UKB_CHIP:/docker/arcas_hla_0.0.1.tar.gz"
        dx_instance_type: "mem1_ssd1_v2_x8"
    }

    output {
        File extracted_regions = "extracted_regions.bam"
        File extracted_regions_index = "extracted_regions.bam.bai"
    }

}
