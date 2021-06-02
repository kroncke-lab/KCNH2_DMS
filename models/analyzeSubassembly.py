#!/usr/bin/python3

# This script inputs a pair of fastq files. It pulls out the barcode from the R2 read,
# then pulls out the full R1 read, then outputs the result in a table.

import sys
import re
import math
inFile1=sys.argv[1]
inFile2=sys.argv[2]
outFile1=sys.argv[3]
mutationPrefix=str(sys.argv[4])
mutationSuffix=str(sys.argv[5])
mutationRef=str(sys.argv[6])
wtRef=str(sys.argv[7])

# Barcode info (R2)
barcodePrefix='GATCT'
barcodeSuffix='CATATGC' 
barcodeStart=29 
barcodeStop=18+barcodeStart

# Mutation info (R1)
# mutationPrefix='AAAGGT'
# mutationSuffix='CTCTGAGTAGC'
# mutationRef= 'TACGGCGCGGCCGTGCTGTTCTTGCTCATGTGC'
mutationStart=len(mutationPrefix)
mutationStop=len(re.split(mutationSuffix, mutationRef)[0])
mutationRevCom=True
mutationRef=re.split(mutationSuffix, mutationRef)[0]
mutationRef=re.split(mutationPrefix, mutationRef)[1]

# Calculate more parameters
barcodePrefixStart=barcodeStart-len(barcodePrefix)
barcodePrefixStop=barcodeStart
barcodeSuffixStart=barcodeStop
barcodeSuffixStop=barcodeStop+len(barcodeSuffix)
mutationPrefixStart=mutationStart-len(mutationPrefix)
mutationPrefixStop=mutationStart
mutationSuffixStart=mutationStop
mutationSuffixStop=mutationSuffixStart+len(mutationSuffix)

       

########### FUNCTIONS ############
def matchesWithNs(refSeq,testSeq):
    for a in range(len(refSeq)):
        refChar=refSeq[a]
        testChar=testSeq[a]
        if refChar!=testChar and testChar!='N':
            return 0
    return 1

def analyzeMutation(ref,mutatedSeq):
    aaDict = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"X", "TAG":"X",
    "TGT":"C", "TGC":"C", "TGA":"X", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}

    firstMutation=True
    mutatedSeq=mutatedSeq[startMutAA:(endAA-startAA)*3+startMutAA]
    if (mutatedSeq in ref):
        return 'wt'
    if 'N' in list(mutatedSeq):
        return 'N'
    for aaNum in range(0,endAA-startAA,1):
        refAA=ref[(startAA+aaNum)*3:(startAA+aaNum)*3+3]
        mutAA=mutatedSeq[aaNum*3:aaNum*3+3]
        if refAA!=mutAA:
            mutation='goodMut_'+aaDict[refAA]+str(aaNum+startAA)+aaDict[mutAA]+'_'+refAA+str(aaNum+startAA)+mutAA
            if firstMutation:
                firstMutation=False
            else:
                return 'Multiple'
    return mutation

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
fh2=open(inFile2,'r')
fh3=open(outFile1,'w')

# Define where the subassembly sequence lines up with the reference
if mutationRevCom:
    mutationRef=revCom(mutationRef)
endAA=(len(wtRef)-len(re.split(mutationRef,wtRef)[1]))//3
startAA=math.ceil(len(re.split(mutationRef,wtRef)[0])/3)
startMutAA=len(re.split(wtRef[startAA*3:3*startAA+6], mutationRef)[0])
 
# Loop through all the lines
keepGoing=1
count=0
barcodePrefixBad=0
barcodeSuffixBad=0
goodBarcodeCount=0
mutationPrefixBad=0
mutationSuffixBad=0
mutationMultipleCount=0
mutationNCount=0
mutationWTCount=0
goodMutationCount=0
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
    b1=fh2.readline().strip()
    b2=fh2.readline().strip()
    b3=fh2.readline().strip()
    b4=fh2.readline().strip()
    
    ###### Analyze barcode #########
    # Make sure barcode has right prefix and suffix    
    if not matchesWithNs(barcodePrefix,b2[barcodePrefixStart:barcodePrefixStop]):
        barcodePrefixBad+=1
        continue
    if not matchesWithNs(barcodeSuffix,b2[barcodeSuffixStart:barcodeSuffixStop]):
        barcodeSuffixBad+=1
        continue
    goodBarcodeCount+=1
    barcode=b2[barcodeStart:barcodeStop]
    a1split=a1.split(' ')
    readName=a1split[0]
    
    ##### Analyze mutation ######
    # Make sure barcode has right prefix and suffix
    mutation=a2[mutationStart:mutationStop]
    if mutationRevCom:
        mutation=revCom(mutation)
        mutationRef=revCom(mutationRef)
    if not (matchesWithNs(mutationPrefix,a2[mutationPrefixStart:mutationPrefixStop])):
        mutationPrefixBad+=1
        mutationType='mutationPrefixBad'
    elif not (matchesWithNs(mutationSuffix,a2[mutationSuffixStart:mutationSuffixStop])):
        mutationSuffixBad+=1
        mutationType='mutationSuffixBad'
    else:
        mutationType=analyzeMutation(wtRef,mutation)

    # Count different types of mutations
    if mutationType[0:4]=='good':
        goodMutationCount+=1
    if mutationType=='Multiple':
        mutationMultipleCount+=1
    if mutationType=='N':
        mutationNCount+=1
    if mutationType=='wt':
        mutationWTCount+=1
    
    # Write output
    output=[readName,barcode,mutation,mutationType,a2]
    outputLine=','.join(output)
    fh3.write(outputLine+'\n')

print(count-1,'reads processed')
print(barcodePrefixBad, 'reads had a bad barcode prefix')
print(barcodeSuffixBad, 'reads had a bad barcode suffix')
print(goodBarcodeCount,'reads had a good barcode and are included in the output file. Here are their details:')
print(mutationPrefixBad, 'reads had a bad mutation prefix')
print(mutationSuffixBad, 'reads had a bad mutation suffix')
print('Of the rest:')
print(goodMutationCount, 'reads had a good single aa mutant')
print(mutationWTCount, 'reads were wt')
print(mutationNCount,'reads had an N in the pilot region')
print(mutationMultipleCount, 'reads had multiple mutations in the pilot region')

fh1.close()
fh2.close()
fh3.close()
