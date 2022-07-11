from Bio.Seq import Seq             #needed for reverse complement
import math                         #needed for logs

controloutput = []
expoutput = []
controlbarcodes = []
count = 0
#input = "test.fq"
controlinput = "CON1_2_combined.fastq"          #filename with control sequencing reads
expinput = "SAR1_2_combined.fastq"              #filenmae with experimental sequenicng reads
controlseqreads = 0
expseqreads = 0

# open control reads file and extract reads, corrected to all have the same strand orientation
with open(controlinput) as fh:

    line = fh.readline()
    while line:
            line = fh.readline().strip()
            count += 1
            if "CCCACTCCTTT" in line:           #checks strand of sequence and makes it reverse complement if in this orientation
                controlseqreads +=1
                line_seq = Seq(line)
                line_seq_rc = line_seq.reverse_complement()
                controloutput.append(str(line_seq_rc))
            if "GAAAGGACGAAACA" in line:
                controlseqreads += 1
                controloutput.append(line)

#search extracted sequence reads for first iteration of CACCG key sequence and extracts following 20 base pairs as the sgRNA sequence
for a in controloutput:
    i = 0
    for element in a:
        i += 1
        if (a[i:i+5]) == "CACCG":
            barcode = a[i+5:i+25]
            #print (barcode)
            controlbarcodes.append(barcode)
            break

# open experimental reads file and extract reads, corrected to all have the same strand orientation
with open(expinput) as fh:

    line = fh.readline()
    while line:
            line = fh.readline().strip()
            count += 1
            if "CCCACTCCTTT" in line:
                expseqreads+=1
                line_seq = Seq(line)
                line_seq_rc = line_seq.reverse_complement()
                expoutput.append(str(line_seq_rc))
            if "GAAAGGACGAAACA" in line:
                expseqreads+=1
                expoutput.append(line)

#search extracted sequence reads for first iteration of CACCG key sequence and extracts following 20 base pairs as the sgRNA sequence

expbarcodes = []

for a in expoutput:
    i = 0
    for element in a:
        i += 1
        if (a[i:i+5]) == "CACCG":
            barcode = a[i+5:i+25]
            expbarcodes.append(barcode)
            break


library = "broadgpp-brunello-library-contents.txt"      #tab delimited text file with sgRNA sequences and other information
summary = {}
readcount = []
# start assembling text for outputfile by adding the column headers
outputfile = "Target Gene ID\tTarget Gene\tTarget Transcript\tGenomic Sequence\t Position of Base After Cut\tStrand\tsgRNA Target Sequence\tTarget Context Sequence\tPAM Sequence\tExon Number\tRule Set 2 score\tControl Reads\tExp Reads\tLogNorm Control\tLogNorm Exp\tLog Fold Change\n"

with open(library) as f:        #opens library file and reads line-by-line
    line = f.readline()
    while line:
        line = f.readline()
        line = line.rstrip()
        linelist = line.split("\t")
        if expbarcodes.count(linelist[6]) != 1111111111111111111111111111111111111111110:       #placeholder - forces program to check all lines. should rewrite.
            explognorm = math.log((expbarcodes.count(linelist[6])/expseqreads*1000000+1),2)     #gets read count for each sgRNA and log2 normalizes it
            controllognorm = math.log((controlbarcodes.count(linelist[6]) / controlseqreads * 1000000 + 1),2)   #log2 normalization of control sgRNAs
            logfoldchange = explognorm-controllognorm                                   #log fold change
            #for each sgRNA, assemble result text with all computations to be appended to end of the line from the initial library file
            result = "\t" + str(controlbarcodes.count(linelist[6])) + '\t' + str(expbarcodes.count(linelist[6]))+"\t"+ str(controllognorm)+"\t" + str(explognorm)+"\t" + str(logfoldchange)
            readcount.append(result)
            outputline = str('\t'.join(linelist) + result + "\n")
            outputfile = outputfile + outputline
        if linelist[6] == "TTTTTCTCACCCGATGAATC":           # catching the end of file to prevent list index error
            break

#appends some details on the number of reads to the end of the file.
outputfile += "\nControl Sequence reads\t"+str(controlseqreads) + "\nExperimental Seq Reads\t" + str(expseqreads)


# write out results file as tab-delimited text

fo = open("results.txt", "w")
fo.write(outputfile)
fo.close()


