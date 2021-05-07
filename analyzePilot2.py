#!/usr/bin/python


# Consider using FastqGeneralIterator
# from Bio.SeqIO.QualityIO import FastqGeneralIterator
'''
>>> from Bio.SeqIO.QualityIO import FastqGeneralIterator
>>> count = 0
>>> total_len = 0
>>> with open("example.fastq") as in_handle:
...     for title, seq, qual in FastqGeneralIterator(in_handle):
...         count += 1
...         total_len += len(seq)
'''

# https://biopython.readthedocs.io/en/latest/chapter_cookbook.html#chapter-cookbook


# This script inputs a fastq file. It pulls out the barcode from the read,
# then pulls out the full R1 read, then outputs the result in a table.

import sys
inFile1=sys.argv[1]
outFile1=sys.argv[2]
# Barcode info (R1)
barcodePrefix='AGATCT' 
barcodeSuffix='CATA'
barcodeStart=20 #python counting
barcodeStop=18+barcodeStart

# Calculate more parameters
barcodePrefixStart=barcodeStart-len(barcodePrefix)
barcodePrefixStop=barcodeStart
barcodeSuffixStart=barcodeStop
barcodeSuffixStop=barcodeStop+len(barcodeSuffix)

########### FUNCTIONS ############
def matchesWithNs(refSeq,testSeq):
    for a in range(len(refSeq)):
        refChar=refSeq[a]
        testChar=testSeq[a]
        if refChar!=testChar and testChar!='N':
            return 0
    return 1

def revCom(seq):
    revCommed=''
    RCdict={'a':'t','c':'g','g':'c','t':'a','A':'T','T':'A','C':'G','G':'C','N':'N','n':'n'}
    for char in seq:
        newChar=RCdict[char]
        revCommed=newChar+revCommed
    return revCommed

########### MAIN ##############

# Open filehandles
fh1=open(inFile1,'r')
fh2=open(outFile1,'w')

# Loop through all the lines
keepGoing=1
count=0
barcodePrefixBad=0
barcodeSuffixBad=0
goodBarcodeCount=0
while keepGoing:
    count+=1
    if count/1000000.0==round(count/1000000.0):
        print(count)
    a1=fh1.readline().strip()
    if len(a1)==0:
        break
    a2=fh1.readline().strip()
    a3=fh1.readline().strip()
    a4=fh1.readline().strip()
    
    ###### Analyze barcode #########
    # Make sure barcode has right prefix and suffix    
    if not matchesWithNs(barcodePrefix,a2[barcodePrefixStart:barcodePrefixStop]):
        barcodePrefixBad+=1
        continue
    if not matchesWithNs(barcodeSuffix,a2[barcodeSuffixStart:barcodeSuffixStop]):
        barcodeSuffixBad+=1
        continue
    goodBarcodeCount+=1
    barcode=a2[barcodeStart:barcodeStop]
    a1split=a1.split(' ')
    readName=a1split[0]
    
    # Write output
    output=[readName,barcode,a2]
    outputLine=','.join(output)
    fh2.write(outputLine+'\n')

percentGood=round(100*goodBarcodeCount/(count-1),1)
print(count-1,'reads processed')
print(barcodePrefixBad, 'reads had a bad barcode prefix')
print(barcodeSuffixBad, 'reads had a bad barcode suffix')
print(goodBarcodeCount,'reads had a good barcode and are included in the output file.')
print(percentGood, '% of reads had a good barcode')

fh1.close()
fh2.close()

