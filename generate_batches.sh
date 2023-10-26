# Batch .cram and .crai files for submission
# https://github.com/DrGBL/HLA_UKB/blob/main/step.02.DNA_Nexus/01.create_batch.sh

# Takes each of 10>60 exome containing files, and batches the 9000 exomes in each into 19 batch files of 500 each. Then combines them back into 19 single files, containing 10k exome files each

#where you want batches to be recorded (on your local machine)
path_local_directory=/well/jknight/users/qmy094/analyses/UKB_CHIP/wdl_pipeline

#where the "Exome OQFE CRAM files" files are located (on your DNAnexus folder)
path_to_cram="UKB_CHIP:/Bulk/Exome sequences/Exome OQFE CRAM files/"

#go to working directory and make the right folders
cd ${path_local_directory}

mkdir -p batches
mkdir -p batches_prelim

#now go to the right location in dna nexus
dx cd ${path_to_cram} 

#list the folders (10 to 60)
dx ls > folders_full_exomes.txt
sed -i 's/\///' folders_full_exomes.txt

## Make batch files (get .cram and .cram.crai for each id's sample)

while read folder; do
dx generate_batch_inputs \
-icram_file='(.*).cram$' \
-icram_file_index='(.*).cram.crai$' \
--path "${path_to_cram}${folder}/" \
-o 'batches_prelim/folder_'${folder}
done < folders_full_exomes.txt

## Then process back into individual files for each batch
while read folder; do
  head -n 1 batches_prelim/folder_${folder}.0000.tsv > batches/batch_${folder}.tsv
  awk 'FNR>1' batches_prelim/folder_${folder}* >> batches/batch_${folder}.tsv
  sed -i 's/cram_file/stage-common.cram_file/g' batches/batch_${folder}.tsv
done < folders_full_exomes.txt


