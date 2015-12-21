import numpy as np
from collections import defaultdict
from ReadsLayoutArranger import ReadsLayoutArranger 
import random

class MatrixComputeEngine:
    refGenome_ = ""
    geneBeginIndList_ = []
    genesCnt_ = 0
    genomeLength_ = 0
    genesLength_ = []
    tumorCnt_ = 0
    genesTumorMaxError_ = np.empty(shape=(genesCnt_, tumorCnt_)) # just for initialization
    genesEarliestReadPositionInFile_ = defaultdict(int)
    MetaDataFetched_ = False
    fileHandler_ = None
    DNA_sequence_filename_ = None
    
    def setFile(self, DNA_sequence_filename):
        self.DNA_sequence_filename_ = DNA_sequence_filename
    
    
    
    def getMetaData(self):
        print '[LOG] function MatrixComputeEngine.getMetaData() invoked'

        self.MetaDataFetched_ = True

        self.fileHandler_ = open(self.DNA_sequence_filename_, 'r')
        
        # 1st line: Genome length
        self.genomeLength_ = int(self.fileHandler_.readline().strip().split()[1])
        
        # 2nd line: Tumor cnt
        self.tumorCnt_ = int(self.fileHandler_.readline().strip().split()[1])
   
        # 3rd line: Reference Genome
        self.refGenome_ = self.fileHandler_.readline().strip().split()[1]
    
        # 4th line: Genes indices (which are separated by space)
        self.geneBeginIndList_ = self.fileHandler_.readline().strip().split(' ')[1:]
        # convert strings to int
        self.geneBeginIndList_ = map(int, self.geneBeginIndList_)
        self.genesCnt_ = len(self.geneBeginIndList_)
        self.genesLength_ = [t - s for s, t in zip(self.geneBeginIndList_, self.geneBeginIndList_[1:])] + [self.genomeLength_ - self.geneBeginIndList_[-1]]
        
        return {'tumor_cnt':self.tumorCnt_, 'genes_cnt': self.genesCnt_}
        
    def computeErrorMatrix(self):
        print '[LOG] function MatrixComputeEngine.computeErrorMatrix() invoked'

        if (not self.MetaDataFetched_):
            self.getMetaData(self.DNA_sequence_filename_)
        
        self.genesTumorMaxError_ =  np.zeros(shape=(self.genesCnt_, self.tumorCnt_), dtype=float)
        
        currentGeneID = 0
        currentGeneEndPosition = self.geneBeginIndList_[1] - 1
        
        self.genesEarliestReadPositionInFile_[currentGeneID] = self.fileHandler_.tell();
        
        readCntList = defaultdict(list)
        mutationCntList = defaultdict(list)
            
        while True:
        #for line in iter(f.readline, ''):
            linePositionInFile = self.fileHandler_.tell()
            line = self.fileHandler_.readline()
            if not line:
                break
            #print '********************************'
            #print 'Line:', line
            # Read
            splits = line.split(' ')
            tumor_id = int(splits[0])
            quality = float(splits[1])
            start = int(splits[2])
            read = splits[3].strip()
            
            if (start > currentGeneEndPosition):
                # we have entered a new gene, so we can ship the information of the old gene
                mut_temp = mutationCntList.pop(currentGeneID)
                read_temp = readCntList.pop(currentGeneID)
                
                # go over the read_temp to make sure there is no zero:
                read_temp = (read_temp == 0) * np.ones(shape=(self.tumorCnt_, self.genesLength_[currentGeneID])) + read_temp
                                
                # size of tumors_max_error_rate is Tx1, where T is the number of tumors
                tumors_max_error_rate = np.amax(mut_temp / read_temp, axis=1)
                self.genesTumorMaxError_[currentGeneID] = tumors_max_error_rate
                
                #print 'For gene', currentGeneID, ': The tumor max error rate is',tumors_max_error_rate
                #print 'And the tumor complete matrix is', mut_temp / read_temp
                print '[Info] Processed gene', currentGeneID
                yield {'gene_ID': currentGeneID, 'tumors_max_error_rate': tumors_max_error_rate}
                currentGeneID += 1
                currentGeneEndPosition = self.geneBeginIndList_[currentGeneID] + self.genesLength_[currentGeneID] - 1
            
            for i, readBase in enumerate(read):
                #print readBase
                positionInGenome = start + i
                #print 'position in genome', positionInGenome
                # First, find which gene this base belongs to
                thisGeneID = -1
                thisGeneEndPosition = -1
                if (positionInGenome <= currentGeneEndPosition):
                    # the base belongs to the current gene
                    thisGeneID = currentGeneID
                    thisGeneEndPosition = currentGeneEndPosition
                else:
                    thisGeneID = currentGeneID + 1
                    thisGeneEndPosition = currentGeneEndPosition + self.genesLength_[thisGeneID]
                    while (positionInGenome > thisGeneEndPosition):
                        thisGeneID += 1
                        thisGeneEndPosition = thisGeneEndPosition + self.genesLength_[thisGeneID]         

                # initialize the data structure if they are empty for that gene
                # (this happens when no read has encountered that gene yet)
                if len(readCntList[thisGeneID])==0:
                    readCntList[thisGeneID] = np.zeros(shape=(self.tumorCnt_, self.genesLength_[thisGeneID]))
                    mutationCntList[thisGeneID] = np.zeros(shape=(self.tumorCnt_, self.genesLength_[thisGeneID]))
                    self.genesEarliestReadPositionInFile_[thisGeneID] = linePositionInFile
                
                positionInGene = positionInGenome - self.geneBeginIndList_[thisGeneID]
                #print 'thisGeneID', thisGeneID, ', tumor_id', tumor_id, ', position in gene', positionInGene
                #print 'size', readCntList[thisGeneID].shape 
                readCntList[thisGeneID][tumor_id, positionInGene] += 1
                if (readBase != self.refGenome_[positionInGenome]):
                    mutationCntList[thisGeneID][tumor_id, positionInGene] += 1
    
        # at this point, there might be some genes still left (the last gene is left, for sure).
        # therefore, these genes should be processed
         
        for geneID in readCntList.keys():
            # go over the read_temp to make sure there is no zero:
            read_temp = (readCntList[geneID] == 0) * np.ones(shape=(self.tumorCnt_, self.genesLength_[geneID])) + readCntList[geneID]
            # size of tumors_max_error_rate is Tx1, where T is the number of tumors
            tumors_max_error_rate = np.amax(mutationCntList[geneID] / read_temp, axis=1)
            self.genesTumorMaxError_[geneID] = tumors_max_error_rate

            
            #print 'For gene', geneID, ': The tumor max error rate is',tumors_max_error_rate
            #print 'And the tumor complete matrix is', mutationCntList[geneID] / read_temp
            print '[Info] Processed gene', geneID
            yield {'gene_ID': geneID, 'tumors_max_error_rate': tumors_max_error_rate}
    
        print '[Info] Done!'
    
    def fetchTumorGeneReads(self, tumorID, geneID):
        print '[LOG] function MatrixComputeEngine.fetchTumorGeneReads(', tumorID, ',', geneID, ') invoked'
        
        filePosition = self.genesEarliestReadPositionInFile_[geneID]
        f = open(self.DNA_sequence_filename_, "r")
        f.seek(filePosition)
        
        geneBeginPos = self.geneBeginIndList_[geneID]
        geneEndPos = geneBeginPos + self.genesLength_[geneID] - 1
        
        while True:
            line = f.readline()
            if not line:
                break
            splits = line.split(' ')
            tumor_id = int(splits[0])
            quality = float(splits[1])
            start = int(splits[2])
            read = splits[3].strip()
            
            real_start_position = start
            
            if (start > geneEndPos):
                break
            if (tumor_id != tumorID):
                continue
            
            readEndPos = start + len(read) - 1
            if (readEndPos < geneBeginPos):
                continue
            
            if (start < geneBeginPos):
                # trim the read if it starts before the gene
                read = read[(geneBeginPos - start):]
                start = geneBeginPos
                
            
            yield {'real_start_position': real_start_position, 'start_position': start, 'read': read, 'quality': quality}
    
    def fetchReferenceGene(self, gene_id):
        print '[LOG] function MatrixComputeEngine.fetchReferenceGene(', gene_id, ') invoked'
        if gene_id > len(self.geneBeginIndList_):
            '[LOG] in fetchReferenceGene: requested gene with id', gene_id, 'does not exist'
            return {'gene': '', 'gene_begin_position': 0}
        gene_pos = self.geneBeginIndList_[gene_id]
        gene_length = self.genesLength_[gene_id]
        
        return {'gene_id': gene_id ,'gene': self.refGenome_[gene_pos : gene_pos + gene_length], 'gene_begin_position': gene_pos}
    
    def getErrorMatrix(self):
        for geneID in range(self.genesCnt_):
            yield {'gene_ID': geneID, 'tumors_max_error_rate': self.genesTumorMaxError_[geneID]}
    
    def getGenesBelowThreshold(self, threshold):
        output = []
        for geneID in range(self.genesCnt_):
            aboveThreshold = False
            for tumorID in range(self.tumorCnt_):
                if self.genesTumorMaxError_[geneID, tumorID] > threshold:
                    aboveThreshold = True
                    break
            if aboveThreshold == False:
                output.append(geneID)
        
        return output
            
        


if __name__ == "__main__":
    DNA_sequence_filename = "../../data/small.data"
    matrix = MatrixComputeEngine(DNA_sequence_filename)
    layout = ReadsLayoutArranger()
    print matrix.getMetaData()
    
    for l in matrix.computeErrorMatrix():
        #print l
        pass
    
    gene_id = 2
    layout.initializeCanvas(matrix.genesLength_[gene_id])
    gene_begin_ind = matrix.geneBeginIndList_[gene_id]    
    for l in matrix.fetchTumorGeneReads(0, gene_id):
        print l
        level = layout.addRead(l['start_position'] - gene_begin_ind, len(l['read']))
        print 'LEVELLLLLLLLL', level
    
    print matrix.fetchReferenceGene(gene_id)
    
