#!/bin/sh

#  Genome_index_script.sh
#  
#
#  Created by Maria Jimenez Ramos on 15/01/2021.
#  


# Downloading reference genome and annotation

cd ~/datasets/genome
wget -O Homo_sapiens_genome.fa.gz ftp://ftp.ensembl.org/pub/release-102/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens_genome.fa.gz
cd ../

cd ~/datasets/annotation
wget -O Homo_sapiens_annotation.gtf.gz ftp://ftp.ensembl.org/pub/release-102/gtf/homo_sapiens//Homo_sapiens.GRCh38.102.gtf.gz
gunzip Homo_sapiens_annotation.gtf.gz
cd ../

#Create the genome indexes

~/star/code/STAR-2.7.7a/bin/Linux_x86_64/STAR --runThreadN 8 \ #number of cores
--runMode genomeGenerate \
--genomeDir ~/datasets/genome_index \
--genomeFastaFiles ~/datasets/genome/Homo_sapiens_genome.fa \
--sjdbGTFfile ~/datasets/annotation/Homo_sapiens_annotation.gtf \
--sjdbOverhang 99 # Determine with QC analysis.
