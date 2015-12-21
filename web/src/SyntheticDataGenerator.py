from random import randint
from random import random
from operator import itemgetter
import os
import numpy as np


bases = ['A', 'C', 'G', 'T']
genomeLength = 2000
tumorCnt = 20
readNum = 2000
readMinLength = 3
readMaxLength = 20
#errorRate = 0.2
minGeneLength = 5
maxGeneLength = 50
outputFileName = 'medium2.data'
healthyGeneMutationRate = 0.08;
tumorGeneMutationRate = 0.2;
healthyGeneRate = 0.7;


script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
abs_file_path = os.path.join(script_dir, outputFileName)


# Generating REFERENCE GENOME 
refGenome = ""
for i in range(0, genomeLength):
    refGenome += bases[randint(0, 3)];
    
print '[Info] Reference genome generated.'    

# specifying the genes
genesStartIndList = []
j = 0
genesStartIndList.append(j)
while (j + maxGeneLength < genomeLength-1):
    geneLength = randint(minGeneLength, maxGeneLength)
    j += geneLength
    genesStartIndList.append(j)

    
# specifying the types of gene (healthy/mutated)
genesCnt = len(genesStartIndList)
genesHealthyState = []
for i in range(genesCnt):
    if (random() < healthyGeneRate):
        # healthy
        genesHealthyState.append(0)
    else:
        #mutated
        genesHealthyState.append(1)


# specifying the mutation rate for each base, based on its healthy state
baseMutationRate = np.zeros(shape=(genomeLength))
geneID = 0
for i in range(genomeLength):
    if (geneID + 1 != len(genesStartIndList) and genesStartIndList[geneID + 1] <= i):
        geneID += 1
        
    if (genesHealthyState[geneID] == 0):
        baseMutationRate[i] = healthyGeneMutationRate
    else:
        baseMutationRate[i] = tumorGeneMutationRate;


print '[Info] Genes indices generated.'    

f = open(outputFileName,'w')

f.write('genome_length ' + str(genomeLength))
f.write('\n')

f.write('tumor_cnt ' + str(tumorCnt))
f.write('\n')

f.write('reference_genome ')
f.write(refGenome)
f.write('\n')

f.write('genes_start_indices ')
f.write(' '.join(str(e) for e in genesStartIndList)) # writing gene start-indices to file
f.write('\n')
print '[Info] Reference genome and genes written to file'     


print '[Info] Now generating reads ... (this may take awhile)'     
# Generating reads
allReads = []
for i in range(0, readNum):
    if (i % 1000 == 0):
        print '[Info] Generating reads', i, 'to', i + 1000
        
    for tID in range(tumorCnt):
        # choose a random length
        readLength = randint(readMinLength, readMaxLength)
        # choose a random starting point
        # make sure that the read doesn't exceed the genome
        start = randint(0, genomeLength - readLength)
        # populate the read
        read = ""
        for j in range(start, start + readLength):
            if (random() < baseMutationRate[j]):
                # add noise to the true base 
                trueBase = refGenome[j]
                ind = bases.index(trueBase)
                # shift the true base 1, 2 or 3 (randomly)
                noisyBase = bases[(ind + randint(1, 3)) % 4]
                read += noisyBase
            else:
                read += refGenome[j]
        
        quality = round(random(),1)
        allReads.append({'tumorID': tID, 'start': start, 'quality': quality,'read': read})

print '[Info] Sorting the reads based on their positions before writing to file (this may take quite awhile)'
sorted_list = sorted(allReads, key=itemgetter('start'))

print '[Info] Writing the reads to file (this may take quite awhile)'

for i, dictEntry in enumerate (sorted_list):
    if (i % 1000000 == 0):
        print '[Info] Writing reads', i, 'to', i + 1000000
        
    f.write(str(dictEntry['tumorID']) + ' ' + str(dictEntry['quality']) + ' ' + str(dictEntry['start']) + ' ' + dictEntry['read'])
    if (i < len(sorted_list) - 1):
        f.write('\n')    
    
f.close() # you can omit in most cases as the destructor will call it

print '[Info] Done!'