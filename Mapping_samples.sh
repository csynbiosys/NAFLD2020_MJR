#!/bin/sh

#  Mapping_samples.sh
#  
#
#  Created by Maria Jimenez Ramos on 15/01/2021.
#  
# Script for mapping the samples to the genome with STAR.

cd ~/dataset/samples
~/star/code/STAR-2.7.7a/bin/Linux_x86_64/STAR # directory to the STAR programme.
for i in {67...75} #Range of numbers in the title of the fastq files.
do
    --runMode alignReads \ #Default mode is alignReads, this is not mandatory.
    --runThreadN 16 \ #number of cores
    --genomeDir ~/datasets/genome_index \
    --readFilesIn ~/datasets/samples/SRR72454${i}_1.fq.gz \ ~/datasets/dataset1/samples/SRR72454${i}_2.fq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ~/dataset/alignment/SRR72454${i} \
    --quantMode GeneCounts #Counts matrix for the alignment file.
done
