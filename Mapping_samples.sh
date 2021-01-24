#!/bin/sh

#  Mapping_samples.sh
#  
#
#  Created by Maria Jimenez Ramos on 15/01/2021.
#  
# STAR will be used for mapping.
cd ~/datasets/#dataset specific for that alignments
~/star/code/STAR-2.7.7a/bin/Linux_x86_64/STAR --genomeLoad LoadAndKeep \
--genomeDir ~/datasets/dataset1/genome_index \
for i in {67...75}
do
    --runMode alignReads \
    --runThreadN 8 \ #number of cores
    --genomeDir ~/datasets/dataset1/genome_index \
    --readFilesIn ~/datasets/dataset1/samples/SRR72454$i_1.fq.gz \ ~/datasets/dataset1/samples/SRR72454_$i_2.fq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate Unsorted \
    --quantMode GeneCounts \
done
--genomeLoad Remove --genomeDir ~/datasets/dataset1/genome_index
