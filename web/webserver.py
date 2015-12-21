import time # for time.sleep
from string import Template
import cherrypy
from src.MatrixComputeEngine import MatrixComputeEngine
from src.ReadsLayoutArranger import ReadsLayoutArranger
from random import randint

import json
import simplejson

from threading import Thread

#fileName = '../data/large.data'

__author__ = 'Dan McDougall <YouKnowWho@YouKnowWhat.com>'

# Trying to cut down on long lines...
jquery_url = 'http://ajax.googleapis.com/ajax/libs/jquery/1.3.2/jquery.min.js'
jquery_ui_url = 'http://ajax.googleapis.com/ajax/libs/jqueryui/1.7.2/jquery-ui.min.js'
jquery_ui_css_url = \
'http://ajax.googleapis.com/ajax/libs/jqueryui/1.7.2/themes/black-tie/jquery-ui.css'

class Webserver(object):
    threshold = 0.5
    fixedGenes = {}
    fileName = None
    computeEngine = MatrixComputeEngine()
    readsLayoutArranger = ReadsLayoutArranger()
    receivedGenes = []


    def __init__(self):
        self.colors = {1: '#0404B4', 2: '#A5DF00', 3: '#8A0808', 4: '#DF01D7'}
        
        

    """An example of using CherryPy for Comet-style asynchronous communication"""
    @cherrypy.expose
    def index(self):
        """Return a basic HTML page with a ping form, a kill form, and an iframe"""
        # Note: Dollar signs in string.Template are escaped by using two ($$)

        html = """\
        
        <!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
        <html xmlns="http://www.w3.org/1999/xhtml" lang="en">
        <head>
            <link rel="stylesheet" type="text/css" href="${jquery_ui_css_url}" media="screen" />
            <link rel="stylesheet" type="text/css" href="css/style.css">
            <script src="https://code.jquery.com/jquery-1.9.1.min.js"></script>
            <script src="javascript/scripts.js"></script>
            <script src="javascript/jquery-ui.js"></script>
            <script src="javascript/jquery-ui.min.js"></script>
            
            
        </head>
        <body>
            <div id="whole_page">
                <div id="main_content">
                    <div id="header">
                        <div id="header_buttons">
                            <form action="/showMatrix" method="post" target="console_iframe" enctype="multipart/form-data">
                                <input type="file" class="headerButton" name="myFile" />
                                <input type="submit" class="headerButton"/>
                            </form>
                        </div>
                        <div id="save_button" style="display: inline-block; width: 30px;" class="headerButton">Save</div>
                    </div>
                    
                    <div id="thresholdDiv" style="display: none;">
                        <form method="post" target="console_iframe" action="/changeThreshold">
                            <label for="threshold_number" style="display: inline-block;">Threshold:</label>
                            <div id = "threshold_number" style="width: 40px; height: 30px; display: inline-block;">0.5</div>
                            <input type="range" name="threshold" id="threshold" value="50" min="0" max="100" data-show-value="true">
                            <input class="submitButton" type="submit" data-inline="true" value="Change">
                        </form>
                    
                        <div class="submitButton" id="hide_gene" onclick="hideGenes();">Hide Genes</div>
                        <div class="submitButton" id="show_gene" onclick="showGenes();">Show Genes</div>
                    </div>
                    
                        
                    <div id = "matrix"></div>
                    <div id = "correlationBox" style="display: none"></div>
                    
                    <div id = "reads_canvas"></div>
                            
                    <iframe name="console_iframe" style="display:none"/> 
                    </div>
            </div>
        </body>
        </html>
        """
        t = Template(html)
        page = t.substitute(
            jquery_ui_css_url=jquery_ui_css_url,
            jquery_url=jquery_url,
            jquery_ui_url=jquery_ui_url)
        return page
    
    
    @cherrypy.expose
    def changeThreshold(self, threshold):
        cherrypy.response.headers['Content-Type'] = 'application/json'
        print '[LOG] Function changeThreshold invoked',
        threshold = float(threshold)
        print 'with threshold =' , threshold
        self.threshold = round(threshold/100,2)
        print self.threshold
        
        for gene_max_error_rate in self.computeEngine.getErrorMatrix():
            print '[LOG] Received from compute engine the mutation array for gene', gene_max_error_rate['gene_ID']
            geneID = gene_max_error_rate['gene_ID']
            if (geneID in self.fixedGenes and self.fixedGenes[geneID] == 1):
                # this gene has already been marked by the user as normal, so do not show it
                continue;
            
            tumors_max_error_rate = gene_max_error_rate['tumors_max_error_rate']
            #tumors_binary_error = tumors_max_error_rate > self.threshold
            list_in_json = json.dumps(tumors_max_error_rate.tolist())
            
            function_invoke_str = '<script>parent.updateCell(' + str(geneID) + ',' + list_in_json + ',' + str(self.threshold) + ')</script>'
            print '[LOG] Invoking JS.parent.updateCell for gene', geneID
            yield function_invoke_str
            
    changeThreshold._cp_config = {'response.stream': True}
    

      
    @cherrypy.expose
    def showMatrix(self, myFile):
        """Incrementally compute the Gene-Tumor matrix"""
        self.computeEngine = MatrixComputeEngine()
        self.readsLayoutArranger = ReadsLayoutArranger()
        self.threshold = 0.5
        self.fixedGenes = {}
        self.fileName = None
        function_invoke_str = '<script>parent.resetDisplay()</script>'
        yield function_invoke_str


        print myFile.filename
        print myFile.content_type
        
        self.fileName = '../data/' + myFile.filename
        self.computeEngine.setFile(self.fileName)
        print '[LOG] Function showMatrix invoked'
        metaData = self.computeEngine.getMetaData()
        print '[LOG] Received from compute engine the medatadata'
              
        genesCnt = metaData['genes_cnt']
        tumorsCnt = metaData['tumor_cnt']
        function_invoke_str = '<script>parent.createMatrix(' + str(genesCnt) + ',' + str(tumorsCnt) + ')</script>'
        print '[LOG] Invoking JS: ', function_invoke_str
        yield function_invoke_str
          
        time.sleep(0.1)
          
        for gene_max_error_rate in self.computeEngine.computeErrorMatrix():
            print '[LOG] Received from compute engine the mutation array for gene', gene_max_error_rate['gene_ID']
            geneID = gene_max_error_rate['gene_ID']
            tumors_max_error_rate = gene_max_error_rate['tumors_max_error_rate']
            self.receivedGenes.append(geneID);
            #tumors_binary_error = tumors_max_error_rate > self.threshold
            list_in_json = json.dumps(tumors_max_error_rate.tolist())
  
            function_invoke_str = '<script>parent.updateCell(' + str(geneID) + ',' + list_in_json + ',' + str(self.threshold) + ')</script>'
            
            print '[LOG] Invoking JS.parent.updateCell for gene', geneID
            yield function_invoke_str
            #time.sleep(1)
    
    # Enable streaming for the ping method.  Without this it won't work.
    showMatrix._cp_config = {'response.stream': True}


    @cherrypy.expose
    def kill_proc(self, **kw):
        """Kill the process """
        return "<strong>Success:</strong> The process was stopped successfully."
    
    @cherrypy.expose
    def fixGeneMutation(self, gene, tumor, newstate):
        cherrypy.response.headers['Content-Type'] = 'application/json'
        print '[LOG] Function fixGeneMutation invoked',
        geneID = int(gene)
        newState = int(newstate)
        print 'with input arguments: gene =', geneID , ', state =', newState        
        self.fixedGenes[geneID] = newState
        return simplejson.dumps({'ret': 0})
    
    @cherrypy.expose
    def getCorrelatedGenes(self, gene):
        cherrypy.response.headers['Content-Type'] = 'application/json'
        print '[LOG] Function getCorrelatedGenes invoked',
        geneID = int(gene)
        print 'with input argument gene =', geneID
        
        if len(self.receivedGenes) > 20:
            time.sleep(1)
            if len(self.receivedGenes) > 100:
                time.sleep(1)
                if len(self.receivedGenes) > 1000:
                    time.sleep(3)
                    
            
        
        relatedGenesIDs = []
        print 'receivedGenes' , self.receivedGenes
        for i in range(3):
            ind = randint(0, len(self.receivedGenes) - 1)
            print 'ind', ind
            gid = self.receivedGenes[ind]
            if gid == geneID:
                continue
            relatedGenesIDs.append(gid)
        
        return simplejson.dumps(relatedGenesIDs)        
            
    
    @cherrypy.expose
    def hideGenes(self):
        cherrypy.response.headers['Content-Type'] = 'application/json'
        genesArray = self.computeEngine.getGenesBelowThreshold(self.threshold);
        return simplejson.dumps(genesArray)
    
    def showGenes(self):
        cherrypy.response.headers['Content-Type'] = 'application/json'
        print '[LOG] Function showGenes invoked'
        
        for gene_max_error_rate in self.computeEngine.getErrorMatrix():
            print '[LOG] Received from compute engine the mutation array for gene', gene_max_error_rate['gene_ID']
            geneID = gene_max_error_rate['gene_ID']
            if (geneID in self.fixedGenes and self.fixedGenes[geneID] == 1):
                # this gene has already been marked by the user as normal, so do not show it
                continue;
            
            tumors_max_error_rate = gene_max_error_rate['tumors_max_error_rate']
            #tumors_binary_error = tumors_max_error_rate > self.threshold
            list_in_json = json.dumps(tumors_max_error_rate.tolist())
            
            function_invoke_str = '<script>parent.updateCell(' + str(geneID) + ',' + list_in_json + ',' + str(self.threshold) + ')</script>'
            
            print '[LOG] Invoking JS.parent.updateCell for gene', geneID
            yield function_invoke_str
            
    changeThreshold._cp_config = {'response.stream': True}

        
    
    @cherrypy.expose
    def showReads(self, gene, tumor):
        cherrypy.response.headers['Content-Type'] = 'application/json'
        print '[LOG] Function showReads invoked',
        geneID = int(gene)
        tumorID = int(tumor)
        print 'with input arguments: gene =' , geneID , ', tumor =', tumorID
        referenceGeneDict = self.computeEngine.fetchReferenceGene(geneID)
        referenceGene = referenceGeneDict['gene']
        geneBeginPosition = referenceGeneDict['gene_begin_position']
        
        self.readsLayoutArranger.initializeCanvas(len(referenceGene))
        
        readsList = []
        for read in self.computeEngine.fetchTumorGeneReads(tumorID, geneID):
            readRelativePosition = read['start_position'] - geneBeginPosition 
            level = self.readsLayoutArranger.addRead(readRelativePosition, len(read['read']))
            read['level'] = level
            readsList.append(read)
        
        #print readsList
        json_obj = json.dumps({'reference_gene': referenceGeneDict , 'reads_list': readsList, 'threshold': self.threshold, 'tumor_id': tumorID})
        return json_obj
#         print '[LOG] Invoking JS.parent.showReadsForGene for gene', geneID
#         function_invoke_str = '<script>showReadsForGene(' + json_obj + ')</script>'
#         yield function_invoke_str
#         
    # Enable streaming for the ping method.  Without this it won't work.
    showReads._cp_config = {'response.stream': True}
    
    
    @cherrypy.expose
    def submit(self, name):
        cherrypy.response.headers['Content-Type'] = 'application/json'
        return simplejson.dumps(dict(title="Hello, %s" % name))
        

cherrypy.config.update({
    'log.screen':True,
    'tools.sessions.on': True,
    #'tools.sessions.locking': 'explicit', # This is the magic bit.
    'checker.on':False,
    'server.socket_host':'localhost'})

#cherrypy.quickstart(StringGenerator(), '/', 'simple.config')
cherrypy.tree.mount(Webserver(), config='webserver.config')
cherrypy.engine.start()
cherrypy.engine.block()