#the following java command only needs to run once, but it needs to be re-ran everytime you change the workflow wdl
#it uses the dna-nexus wdl compiling java application, available here: https://github.com/dnanexus/dxCompiler/releases
#note the name of my project here was "exome_full", change this to whatever is needed for you
#it will create a workflow in the dna-nexus directory "my_workflow"

#now move to the local directory where you'll work (in my case, where the "batches" folder is located, created from step 01)
path_local_directory=/well/jknight/users/qmy094/analyses/UKB_CHIP/wdl_pipeline
cd ${path_local_directory}

java -jar dxCompiler-2.5.0.jar compile CHIP_processing.wdl  -project UKB_CHIP -folder CHIP_processing/workflow/

#on dna nexus, just go at the root of your project
dx cd '/'

#now you're ready to call the workflow
#here I hive an example where I run batch 10. It needs to be done for 10 to 60 in order to process the entire UKB

batch=10

dx run CHIP_processing/workflow/CHIP_processing \           
  --batch-tsv batches/batch_${batch}.tsv \    #the batch, this is the result of 01.create_batch.sh
  -istage-common.CHIP_regions="UKB_CHIP:CHIP_references/CHIP_genes_38.bed" \
  -istage-common.ref_genome="UKB_CHIP:CHIP_references/38/Homo_sapiens_assembly38.fasta.gz" \      #the next four files need to be on your project on dna-nexus, i.e. NOT on your local cluster
  -istage-common.ref_genome_index="UKB_CHIP:CHIP_references/38/Homo_sapiens_assembly38.fasta.gz.fai" \
  -istage-common.ref_genome_dict="UKB_CHIP:CHIP_references/38/Homo_sapiens_assembly38.dict" \
  -istage-common.ref_genome_gzi="UKB_CHIP:CHIP_references/38/Homo_sapiens_assembly38.fasta.gz.gzi" \
  --priority low \            #priority for the workers on dna-nexus, see here: https://dnanexus.gitbook.io/uk-biobank-rap/working-on-the-research-analysis-platform/managing-job-priority
  --batch-folders \           #will create a different folder for every result in your batch
  --destination=UKB_CHIP:CHIP_processing/output/batch_10_CHIPextract/output_${batch}      #the destination of the results (you do not need to create it yourself, dna-nexus will create this folder automatically)
