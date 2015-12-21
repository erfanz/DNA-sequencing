import numpy as np

class ReadsLayoutArranger:
    canvasMatrix = None
    geneLength = 0
    maxLayoutLevels = 1000
    
    def initializeCanvas(self, gene_length):
        self.geneLength = gene_length
        self.canvasMatrix = np.ndarray(shape=(self.geneLength, self.maxLayoutLevels), dtype=bool)
        self.canvasMatrix.fill(False)
        
    def addRead(self, relative_start_position, read_length):
        if (relative_start_position < 0):
            # trim the read
            read_length = read_length + relative_start_position
            relative_start_position = 0
        
        #print 'relative', relative_start_position, ', readlength', read_length
        # first, find the first vertical position of read with respect to already laid out reads    
        for level in range(self.maxLayoutLevels):
            #print 'Level: ', level
            if self.doesFit(relative_start_position, level):
                # this place is empty. the read will be placed here
                for i in range(relative_start_position, relative_start_position + read_length):
                    #print 'i', i, ', gene length', self.geneLength , ', relative start', relative_start_position
                    
                    if i >= self.geneLength:
                        break
                    self.canvasMatrix[i, level] = True
                
                return level
        
        # all the levels are full. A new level must be added
        # TODO: needs to be fixed
        return self.maxLayoutLevels - 1
    
    def doesFit(self, position, level):
        if (not self.canvasMatrix[position, level]) \
            and (position == 0 or (position > 0 and not self.canvasMatrix[position - 1, level])):
            return True
        else:
            return False
         
    
        
