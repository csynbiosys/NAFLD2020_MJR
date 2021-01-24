#!/bin/sh

#  Generating_signal.sh
#  
#
#  Created by Maria Jimenez Ramos on 21/01/2021.
#  
#Generate the signal file output to visualize the alignments with programs like IGV or UCSC genomic browser.

~/star/code/STAR-2.7.7a/bin/Linux_x86_64/STAR --runMode inputAlignmentsFromBAM \
--inputBAMfile Aligned.sortedByCoord.out.bam \ #name of the file
--outWigType bedGraph \
--outWigStrand Unstranded

