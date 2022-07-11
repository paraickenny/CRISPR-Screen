# CRISPR-Screen
Tools for analysis of next-generation sequence data from CRISPR screens

Tool developed for the Human CRISPR Knockout Pooled Library (Brunello), available from AddGene.
https://www.addgene.org/pooled-library/broadgpp-human-knockout-brunello/

Files:

1. broadgpp-brunello-library-contents.txt
Tab-delimited text file from Addgene containing data on all lentiviral constructs in the library, including sgRNA insert and the associated gene.

2. FASTQ-to-results.py
This imports FASTQ sequence data from two files (presumably a 'control' and 'experimental' or some other relevant comparison).
Each sequence is checked to determine it's orientation. Sequences opposite the expected orientation are reverse-complemented.
Based on the key sequence, the sgRNA insert sequence for each FASTQ sequence is extracted.
For each sgRNA in the broadgpp-brunello-library-contents.txt file, the number of instances of that sgRNA in the control and experimental files is determined.
Raw sequence counts are normalized to the number of reads for each library, log2 transformed and the fold-difference between control and experimental for each sgRNA is determined. 
For reach sgRNA, these counts and calculated values are appended to the appopriate line in the broadgpp-brunello-library-contents.txt file and the combined data are written out as results.txt
