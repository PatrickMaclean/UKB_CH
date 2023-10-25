version 1.0

#need to change docker image in the hla_calling task runtime section

workflow hla_calling_wf {
    input {
        File cram_file
        File cram_file_index 
        File CHIP_regions
        File ref_genome         #https://github.com/broadinstitute/gatk/raw/master/src/test/resources/large/Homo_sapiens_assembly38.fasta.gz
        File ref_genome_index   #https://github.com/broadinstitute/gatk/raw/master/src/test/resources/large/Homo_sapiens_assembly38.fasta.gz.fai
        File ref_genome_dict    #https://github.com/broadinstitute/gatk/raw/master/src/test/resources/large/Homo_sapiens_assembly38.dict
        File ref_genome_gzi     #https://github.com/broadinstitute/gatk/raw/master/src/test/resources/large/Homo_sapiens_assembly38.fasta.gz.gzi
    }

    call cram_extract { 
        input: cram_file = cram_file, cram_file_index = cram_file_index, CHIP_regions = CHIP_regions, ref_genome = ref_genome, ref_genome_index = ref_genome_index, ref_genome_dict = ref_genome_dict, ref_genome_gzi = ref_genome_gzi
    }

    call bam_index {
        input: CHIP_regions = cram_extract.CHIP_regions
    }

    output {
        File CHIP_regions = cram_extract.CHIP_regions
        File CHIP_regions_index = bam_index.CHIP_regions_index
    }
  
}

task cram_extract {
    
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
        time samtools view --reference ~{ref_genome} -L ~{CHIP_regions} -o CHIP_regions.bam ~{cram_file} 

    >>>

    runtime {
        docker: "dx://UKB_CHIP:/docker/arcas_hla_0.0.1.tar.gz"
        dx_instance_type: "mem1_ssd1_v2_x8"
    }

    output {
        File CHIP_regions = "CHIP_regions.bam"
    }

}

task bam_index {
    input {
        File CHIP_regions
    }

    command <<<

        samtools index ~{CHIP_regions} 
    >>>

    runtime {
        docker: "dx://UKB_CHIP:/docker/arcas_hla_0.0.1.tar.gz"
        dx_instance_type: "mem1_ssd1_v2_x8"
    }

    output {
        File CHIP_regions_index = "CHIP_regions.bam.bai"
    }

}
