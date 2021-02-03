#!/bin/sh

#  QC_samples.sh
#  
#
#  Created by Maria Jimenez Ramos on 18/01/2021.
#  

cd ~/datasets/samples #Folder where the fastq files are saved.
#Perform quality control in each sample.
parallel fastqc {} -o ../results ::: *.fastq.gz

