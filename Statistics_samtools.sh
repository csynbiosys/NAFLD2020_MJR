#!/bin/sh

#  Statistics_samtools.sh
#  
#
#  Created by Maria Jimenez Ramos on 21/01/2021.
#  
#Create the index for the bam file that can be used as input on IGV to visualise the alignments.
samtools index -b SRR7245467Aligned.sortedByCoord.out.bam #it generates a .bai file

#Statistics and analysis about the alignement files generated with STAR.
samtools stats SRR7245467Aligned.sortedByCoord.out.bam #Collects statistics from the BAM files and outputs in a text format. Similar output to the one produced in Log.final.out by STAR, but is more detailed.

samtools flagstat SRR7245467Aligned.sortedByCoord.out.bam #does a full pass through the input file to calculate and print statistics. It gives counts of the FLAG categories with a QC pass and QC fail.
