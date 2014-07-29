'''
Created on 18.07.2014

@author: Axel Wilbertz
'''
import hashlib
import os
from copy import deepcopy

class FileObject(object):
    '''
    classdocs
    '''
    filepath = "";
    genecolumn=0;
    delimiter = "";
    values = [];
    titles = [];
    numrows = 0;
    numcols = 0;


    def __init__(self,filepath,delimiter):
        '''
        Constructor
        '''
        self.filepath = filepath;
        self.delimiter = delimiter;
    
    def setFilepath(self,filepath):
        self.filepath = filepath;
        
    def setGeneColumn(self,id):
        self.genecolumn = id;
        
    def getGeneColumn(self):
        return self.genecolumn;
    
    
    def findGeneColumn(self,path):
        f = open('maps/GeneColumn.tsv','r')
        d = {}
        for l in f:
            if l != None or l != '':
                line = l.strip().split('\t')        
                d.setdefault(line[0],line[1])
        f.close()
        if d.setdefault(os.path.basename(self.filepath)):
            return d.setdefault(os.path.basename(self.filepath))
        else:
            raise ValueError('no gene column found check maps/GeneColumn.tsv and enter filename and gene column')
    
    def getGenePos(self,columnname):
        if self.titles:
            self.genecolumn = self.titles.index(columnname)
            return self.genecolumn
        else:
            f = open(self.filepath,'r')
            title = f.readline().strip().split(self.delimiter)
            self.genecolumn = title.index(columnname)
            return self.genecolumn
            
    
    def getTitles(self):
        f = open(self.filepath,'r')
        self.titles = f.readline().strip().split(self.delimiter)
        f.close()
        return self.titles
    
    def getOffsetValues(self,offset,gene):
        f = open(self.filepath,'r')
        f.seek(offset)
        i = 0
        match = False        
        #this is ugly better use an temparray vor appending and set values to temparray at the end
        print self.genecolumn

        self.values=None        
        self.values=[]
        
        for line in f:
            line = line.split('\t')
            if line[self.genecolumn] == gene:
                match = True
                self.values.append(line)
                
            if self.values and line[self.genecolumn] != gene and i>500:
                print 'Stopped search Gene was found'
                f.close()
                break
            if match:
                i=i+1  

        f.close()
        return self.values
    
    def getAllValues(self):
        tempval = []
        f = open(self.filepath,'r')
        f.readline()
        for line in f:
            line = line.strip().split(self.delimiter)
            tempval.append(line)
         
        f.close()    
        self.values = tempval        
        return tempval
    
    def getGeneValues(self,gene,multi):
        tempval = [];
        f = open(self.filepath,'r');
        for line in f:
            line = line.split(self.delimiter)
            if line[self.genecolumn]==gene and not multi:
                # for correlation plot only one line should be found. for splicing plot there can be several transcripts for one gene
                tempval=line
            elif line[self.genecolumn]==gene and multi:
                tempval.append(line)
        f.close()
        if not tempval:
            raise ValueError('No Values found. Check Gene Col Position')
        else:                
            self.values = tempval; 
            return tempval;               
        
    def countRowsNCols(self,filepath):
        f = open(filepath, 'r')
        col = f.readline().count('\t')
        row = 1
        for line in f:
            row += 1        
        f.close()         
        # da 5 messwerte vorhanden sind allerdings nur 4 Tabstopps dazwischen.    
        self.numcols=col+1
        self.numrows= row
        return [col + 1, row]

    
    def printFilepathLinewise(self):
        f = open(self.filepath,'r')
        for line in f:
            print line
        f.close()



class Mapper(FileObject):
    '''
    Maps the Patients(colum titles) to the related MDS Subtype which can be found in the KorrelationMapping File
    '''
    mapperdict = {}
    positions = [[],[]]
    def __init__(self,filepath,delimiter):
        '''
        Constructor
        '''
        self.filepath = filepath;
        self.delimiter = delimiter;
    
    def extractGroupMapping(self):
        f = open(self.filepath, 'r')
        dict = {}
        for line in f:
            line = line.strip().split(self.delimiter)
            dict.setdefault(line[1], []).append(line[0])
    
        f.close()
        self.mapperdict = dict
        return dict
    
    def findPositions(self,titles,mapperdict):
        pos = [[],[]]

        for k in mapperdict:
            pos[0].append(k)
            temp = []
            for v in mapperdict.get(k):
                temp.append(titles.index(v))
            pos[1].append(temp)
        self.positions = pos
        return pos
    
class RefGene(FileObject):

    indexfilepath = ""
    transpos = -1
    refstartpos = -1
    refstoppos = -1
    expstartpos = -1
    expstoppos = -1
    
    refgenedict = {}
    refgenestartstoplist=[]
    genelist = []
    
    exon =[]
    intron = []
    
    def __init__(self,filepath,indexfilepath,genecolumn,delimiter,transpos,refstartpos,refstoppos,expstartpos,expstoppos,titles, values):
        self.filepath=filepath
        self.indexfilepath=indexfilepath
        self.genecolumn=genecolumn
        self.delimiter=delimiter
        self.transpos=transpos
        self.refstartpos=refstartpos
        self.refstoppos=refstoppos
        self.expstartpos=expstartpos
        self.expstoppos=expstoppos
        self.values=values
        
    def extractRefGeneIndexFileEntries(self):
        f = open(self.indexfilepath, 'r')
        f.readline()
        d = {}
        i = 0
        for line in f:
            l = line.strip().split('\t')
            d.setdefault(l[0], l[1])
            i = i + 1
        f.close()
        self.refgenedict=d
        return d
    
    def extractRefGeneValues(self):
        
        refgene = open(self.filepath, 'r')
        genelist=[]
        #ersten zwei anfangsbuchstaben des gens der genexprdatei in der indexdatei suchen und dorthin springen
        print self.values[0][self.genecolumn][:2]
        off = self.refgenedict.setdefault(self.values[0][self.genecolumn][:2])
        refgene.seek(long(off))
        print 'offset ' + str(off)
        #solange die ersten beiden buchstaben aus genexpr und refgene gleich sind geht er die refgene durch  
        refline = refgene.readline().strip().split('\t')
        while self.values[0][self.genecolumn][:2] == refline[12][:2] :
            if self.values[0][self.genecolumn] == refline[12]:            
                #print 'expr ' + values[0][genepos] + ' ref ' + refline[12] 
                genelist.append(refline)
            elif self.values[0][self.genecolumn] != refline[12] and genelist:
                break
            refline = refgene.readline().strip().split('\t')  
              
        else:
            if not genelist:
                print 'Gene ' + self.values[0][self.genecolumn] + ' not found in refGene file'
                refgene.close()
                return None   
        refgene.close()    
                    
        #zieht sich aus der probesets den start und stop bereich
        refgenestartstoplist = []
        for v in genelist:
            refstart = v[self.refstartpos].split(',')
            refstop = v[self.refstoppos].split(',')
            trans = v[self.transpos]
            refstart = filter(None, refstart)
            refstop = filter(None, refstop)                   
            refgenestartstoplist.append([refstart, refstop,trans])
        #geht die transkript liste durch ueberprueft jedes probeset ob es in dem bereich liegt
        self.refgenestartstoplist=refgenestartstoplist
        self.genelist=genelist
        return refgenestartstoplist
        
        
    def calcIntronExon(self):
        #checkIndexFile('test/refGene.tsv')
        #LOC256021 multiple entries in refgene
           
        exon = []
        intron = []
        
        for t in self.refgenestartstoplist:           
            tempe = []
            tempi=[]
            j=0
            for p in self.values:
                #print 'r '+str(p[expstartpos])+' '+str(p[expstoppos])
                i = 0            
                while i < len(t[0]):
                    if long(t[0][i]) <= long(p[self.expstartpos]) and  long(p[self.expstoppos]) <= long(t[1][i]):                     
                        #print 'exon refstart ' +str(t[0][i]) +' probestart '+ str(p[14])+' probestop '+str(p[15])  +' refstop '+str(t[1][i])            
                        tempe.append([j,t[0][i]])
                    
                    elif i < len(t[0])-1: 
                        #print 'intron refstart ' +str(t[1][i]) +' probestart '+ str(p[14])+' probestop '+str(p[15])  +' refstop '+str(t[0][i+1])            
                        if long(t[1][i]) <= long(p[self.expstartpos]) and long(p[self.expstoppos]) <=  long(t[0][i+1]):                        
                            tempi.append([j,t[1][i]])               
                    i=i+1
                j=j+1
            exon.append(tempe)
            intron.append(tempi)
            
        self.exon = exon
        self.intron = intron
        print '\nintrons'
        print self.intron
        print 'exons'
        print self.exon
    
        return [exon,intron]
    
class IndexFile(FileObject):
    
    hashvalue = -1
    indexfilepath =''
    filepath=''
    
    def __init__(self,filepath,genecolumn):
        self.filepath=filepath
        self.genecolumn = genecolumn
        
   
    def checkIndexFile(self):
        """
        calculates the indexfilepath. checks if an indexfile already exists. creates new indexfile if none exists. compares indexfiles hashvalues
        """
        print 'checkIndexFile filepath '+str(self.filepath)    
        
        self.indexfilepath = self.pathToindexpath()
        
        print 'checkIndexFile indexfilepath '+str(self.indexfilepath)
        
        try:
            f = open(self.indexfilepath,'r')
            md5check = self.md5Checksum()
            md5file = f.readline().strip()
            #print 'mdcheck '+str(mdcheck)+' mdfile '+str(mdfile)
            if md5check == md5file:
                print 'file '+str(self.indexfilepath)+ ' exists and hashes equal.'            
                f.close()
                return True
            else:
                print 'file '+str(self.indexfilepath)+ ' exists but wrong hash value.'
                raise IOError
            
        except IOError:
            print 'file '+str(self.indexfilepath)+ ' does not exist. creating file'
            #self.createIndexFile(self.filepath,self.genecolumn,self.delimiter)
            self.createIndexFile()
            return True

    
    def md5Checksum(self):
        with open(self.filepath, 'rb') as fh:
            m = hashlib.md5()
            while True:
                data = fh.read(8192)
                if not data:
                    break
                m.update(data)
            return m.hexdigest()

    def getIndexOffset(self,gene):
        if not gene or not self.filepath:
            return None
        else:
                
            f = open(self.indexfilepath,'r')
            f.readline()
            for l in f:
                line = l.strip().split('\t')
                if gene[:2] == line[0]:
                    f.close()
                    return line[1]
            f.close()
            return None    

    def createIndexFile(self):
        print 'filepath '+str(self.filepath)
        print 'genepos '+str(self.genecolumn) 
        
        f = open(self.filepath,'r')
        #title = self.titles
        preoffset = f.tell()
        preline=f.readline().split('\t')
        postline=deepcopy(preline)
        postoffset = f.tell()
        temp = []
        index = []
        try:
            temp.append(preline[self.genecolumn][:2])
        except TypeError:
            print 'Please check if file is entered in the maps/GeneColumn.tsv file. Otherwise Gencolumn can not be found'
            return None
        temp.append(str(preoffset)+'\n')
        index.append(temp)
        #calculate index
        try:
            while True:
                preoffset = f.tell()
                line = f.readline()
                if(line == ''):
                    break
                preline = line.split('\t')
    
                
                if preline[self.genecolumn][:2] != postline[self.genecolumn][:2]:
                    temp = []
                    temp.append(preline[self.genecolumn][:2])
                    temp.append(str(preoffset)+'\n')
                    print temp
                    index.append(temp)                                   
                
                postline = deepcopy(preline)
                postoffset = long(preoffset)
    
        except IndexError:
                'Index Error'
    
        f.close()
        print 'index und hashwert berechnen'
        #write hash value in first line of file
        index.insert(0,str(self.md5Checksum())+'\n')    
        #indexfilepath = self.pathToindexpath()
        
        f = open(self.indexfilepath,'w')
        f.write(index[0])
        for v in index[1:]:
            f.write('\t'.join(v))
        f.close
        print 'file '+self.indexfilepath+' created'
        return self.indexfilepath
    

    def pathToindexpath(self):
        filename = os.path.basename(self.filepath)
        dirname = os.path.dirname(self.filepath)
        indexfilepath = 'indexfiles\\'+filename.split('.')[0]+'_index.tsv'
        self.indexfilepath=indexfilepath
        return indexfilepath

    