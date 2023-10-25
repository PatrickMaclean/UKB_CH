version 1.0
# Source https://dnanexus.gitbook.io/uk-biobank-rap/science-corner/guide-to-analyzing-large-sample-sets#tip-submitting-the-large-batch-of-jobs-via-swiss-army-knife

task CRAM_extract{
    input {
        Array[File]+ cram_file
        File ref_genome
        File CHIP_regions
    }

    command <<<
        set -x -e -o pipefail
        mkdir output
        for input in ~{sep=" " mapped_read}; do
            file_prefix=$( basename $input ".cram")
        time samtools view --reference ~{ref_genome} \
        -L ~{CHIP_regions} \
        -o "~{file_prefix}_CHIPregions.bam" \
        ~{cram_file} && \
        samtools index "~{file_prefix}_CHIPregions.bam"
        done
    >>>

    output {
        Array[File] CHIPregionsextracted = glob("output/*.CHIPregions.bam")
        Array[File] CHIPregionsextractedindex = glob("output/*.CHIPregions.bam.bai")
    }
    runtime {
        docker: "dx://UKB_CHIP:/docker/arcas_hla_0.0.1.tar.gz"
        dx_instance_type: "mem1_ssd1_v2_x8"
    }

    parameter_meta {
    cram_file: {
        description: "mapped short read data",
        patterns: ["*.cram"],
        stream: true
    }
    ref_genome: {
        description: "reference genome",
        patterns: ["*.fa","*.fasta"]
    }
    CHIP_regions: {
        description: "list of CHIP gene regions in bed format",
        patterns: ["*.bed"]
    }
    }
}
