'''
Created on 23.04.2013
@author: naxobIdeaPad
'''
import wx
import os
import wx.grid
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from numpy import *
from copy import deepcopy
import weakref
import math
import hashlib
from copy import deepcopy

def md5Checksum(filePath):
    with open(filePath, 'rb') as fh:
        m = hashlib.md5()
        while True:
            data = fh.read(8192)
            if not data:
                break
            m.update(data)
        return m.hexdigest()

def getGeneColMatch(filepath):
    f = open('maps/GeneColumn.tsv','r')
    d = {}
    for l in f:
        if l != None or l != '':
            line = l.strip().split('\t')        
            d.setdefault(line[0],line[1])
    f.close()
    return d.setdefault(os.path.basename(filepath))

def getExon(genepos, values):
    #checkIndexFile('test/refGene.tsv')
    #LOC256021 multiple entries in refgene
    refgene = open('test/refGene.tsv', 'r')
    f3 = open('indexfiles/refGene_index.tsv', 'r')
    f3.readline()
    d = {}
    i = 0
    for line in f3:
        l = line.strip().split('\t')
        d.setdefault(l[0], l[1])
        i = i + 1
    f3.close()
    genelist = []
    
                     
    #ersten zwei anfangsbuchstaben des gens der genexprdatei in der indexdatei suchen und dorthin springen
    off = d.setdefault(values[0][genepos][:2])
    refgene.seek(long(off))
    print 'offset ' + str(off)
    #solange die ersten beiden buchstaben aus genexpr und refgene gleich sind geht er die refgene durch  
    refline = refgene.readline().strip().split('\t')
    while values[0][genepos][:2] == refline[12][:2] :
        if values[0][genepos] == refline[12]:            
            #print 'expr ' + values[0][genepos] + ' ref ' + refline[12] 
            genelist.append(refline)
        elif values[0][genepos] != refline[12] and genelist:
            break
        refline = refgene.readline().strip().split('\t')  
          
    else:
        if not genelist:
            print 'Gene ' + values[0][genepos] + ' not found in refGene file'
            refgene.close()
            return None   
    refgene.close()    
                
    #zieht sich aus der probesets den start und stop bereich
    refgenestartstoplist = []
    for v in genelist:
        refstart = v[9].split(',')
        refstop = v[10].split(',')
        trans = v[1]
        refstart = filter(None, refstart)
        refstop = filter(None, refstop)                   
        refgenestartstoplist.append([refstart, refstop,trans])
    #geht die transkript liste durch ueberprueft jedes probeset ob es in dem bereich liegt
    exon = []
    intron = []
    for t in refgenestartstoplist:
        tempe = []
        tempi=[]
        j=0
        for p in values:
            i = 0            
            while i < len(t[0]):
                if long(t[0][i]) <= long(p[14]) and  long(p[15]) <= long(t[1][i]):                     
                    #print 'exon refstart ' +str(t[0][i]) +' probestart '+ str(p[14])+' probestop '+str(p[15])  +' refstop '+str(t[1][i])            
                    tempe.append([j,t[0][i]])
                
                elif i < len(t[0])-1: 
                    #print 'intron refstart ' +str(t[1][i]) +' probestart '+ str(p[14])+' probestop '+str(p[15])  +' refstop '+str(t[0][i+1])            
                    if long(t[1][i]) <= long(p[14]) and long(p[15]) <=  long(t[0][i+1]):                        
                        tempi.append([j,t[1][i]])               
                i=i+1
            j=j+1
        exon.append(tempe)
        intron.append(tempi)
    
    return [exon,intron]
        

def checkIndexFile(filepath):
    print 'checkIndexFile filepath '+str(filepath)    
    
    indexfilepath = pathToindexpath(filepath)
    
    print 'checkIndexFile indexfilepath '+str(indexfilepath)
    try:
        f = open(indexfilepath,'r')
        md5check = md5Checksum(filepath)
        md5file = f.readline().strip()
        #print 'mdcheck '+str(mdcheck)+' mdfile '+str(mdfile)
        if md5check == md5file:
            print 'file '+str(indexfilepath)+ ' exists and hashes equal.'            
            f.close()
            return True
        else:
            print 'file '+str(indexfilepath)+ ' exists but wrong hash value.'
            raise IOError
        
    except IOError:
        print 'file '+str(indexfilepath)+ ' does not exist. creating file'
        createIndexFile(filepath,getgenpos(str(filepath),getGeneColMatch(filepath),'\t'))
        return True

def getIndexOffset(gene,filepath):
    if not gene or not filepath:
        return None
    else:
        
        indexfilepath = pathToindexpath(filepath)

        f = open(indexfilepath,'r')
        f.readline()
        for l in f:
            line = l.strip().split('\t')
            if gene[:2] == line[0]:
                f.close()
                return line[1]
        f.close()
        return None    

                   
def absolutePathtoRelative (path):  
    path = str(path)
    path = path.encode('string-escape')
    path = path.split('\\')
    path = path[(len(path)-3)]+'/'+path[len(path)-1]
    return path

def pathToindexpath(path):
    """
    indexfilepath=path.split('.')
    indexfilepath[0]=indexfilepath[0]+'_index'
    indexfilepath='.'.join(indexfilepath)

    indexfilepath=indexfilepath.split('/')
    indexfilepath[0]='indexfiles'
    indexfilepath='/'.join(indexfilepath)
    """
    filename = os.path.basename(path)
    dirname = os.path.dirname(path)
    indexfilepath = dirname+'\\indexfiles\\'+filename.split('.')[0]+'_index.tsv'
    return indexfilepath

def createIndexFile(filepath,genepos):
    print 'filepath '+str(filepath)
    print 'genepos '+str(genepos) 
    
    f = open(filepath,'r')
    head = f.readline()
    preoffset = f.tell()
    preline=f.readline().split('\t')
    postline=deepcopy(preline)
    postoffset = f.tell()
    temp = []
    index = []
    try:
        temp.append(preline[genepos][:2])
    except TypeError:
        print 'Please check if file is entered in the maps/GeneColumn.tsv file. Otherwise Gencolumn can not be found'
        return None
    temp.append(str(preoffset)+'\n')
    index.append(temp)

    try:
        while True:
            preoffset = f.tell()
            line = f.readline()
            if(line == ''):
                break
            preline = line.split('\t')

            
            if preline[genepos][:2] != postline[genepos][:2]:
                temp = []
                temp.append(preline[genepos][:2])
                temp.append(str(preoffset)+'\n')
                print temp
                index.append(temp)
                #print '\npostline '+str(postline)
                #print 'preline '+str(preline)
                #print 'temp '+str(temp)
                ##print 'index '+str(index)
                               
            
            postline = deepcopy(preline)
            postoffset = long(preoffset)

    except IndexError:
            'Index Error'

    f.close()
    print 'index und hashwert berechnen'
    index.insert(0,str(md5Checksum(filepath))+'\n')

    """
    filepath=filepath.split('.')
    filepath[0]=filepath[0]+'_index'
    filepath='.'.join(filepath)
    
    filepath=filepath.split('/')
    filepath[0]='indexfiles'
    filepath='/'.join(filepath)
    """

    indexfilepath = pathToindexpath(filepath)
    f = open(indexfilepath,'w')
    f.write(index[0])
    for v in index[1:]:
        f.write('\t'.join(v))
    f.close
    print 'file '+indexfilepath+' created'


def getgenpos(array, gencolname, *delimiter):      
        if type(array) == str and array.index('.'):
            file = open(array, 'r')            
            array = file.readline()
            array = array.split(delimiter[0])   
            file.close()         
        if(type(gencolname) == str):
            if type(array[0]) == list:
                return array[0].index(gencolname)
            else:
                return array.index(gencolname)
        elif type(gencolname) == list :
            indexarray = []
            if type(array[0]) == list:            
                for elem in gencolname:           
                    indexarray.append(array[0].index(elem))
            else:
                for elem in gencolname:           
                    indexarray.append(array.index(elem))
            return indexarray  

#not in use
def assignIDs(list):
    '''Take a list of strings, and for each unique value assign a number.
    Returns a map for "unique-val"->id.
    '''
    sortedList = sorted(list)

    #taken from
    #http://stackoverflow.com/questions/480214/how-do-you-remove-duplicates-from-a-list-in-python-whilst-preserving-order/480227#480227
    seen = set()
    seen_add = seen.add
    uniqueList = [ x for x in sortedList if x not in seen and not seen_add(x)]

    return  dict(zip(uniqueList, range(len(uniqueList))))

#not in use
def plotData(inData, color):
    x, y = zip(*inData)

    xMap = assignIDs(x)
    xAsInts = [xMap[i] for i in x]


    plt.scatter(xAsInts, y, color=color)
    plt.xticks(xMap.values(), xMap.keys())


def getGroupMapping(filepath):
    f = open(filepath, 'r')
    dict = {}
    for l in f:
        l = l.strip().split(':')
        dict.setdefault(l[1], []).append(l[0])

    f.close()
    return dict

def countRowsNCols(filepath):
        f = open(filepath, 'r')
        col = f.readline().count('\t')
        row = 1
        for line in f:
            row += 1        
        f.close()         
        # da 5 messwerte vorhanden sind allerdings nur 4 Tabstopps dazwischen.     
        return [col + 1, row]

class ChooseFrame(wx.Frame):

    title = "Choose Genes and Colors"    
    def __init__(self, parent, values, *args, **kwargs):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, title=self.title)
        self.panel = DemoPanel(self)
        self.grid = wx.grid.Grid(self.panel)                
        self.grid.CreateGrid(len(values[0]), 2)
        print values
        mdict = getGroupMapping("maps/groupmapping.csv")
        values.append(mdict)
        i = 0
        for v in values[0]:       
            if mdict.has_key(v):
                self.grid.SetCellValue(i, 1, mdict.get(v))                     
            self.grid.SetRowLabelValue(i, v)
            self.grid.SetCellValue(i, 0, values[1][i])
            i = i + 1         
        
        print mdict
        self.panel.box.Add(self.grid, 2, wx.EXPAND)
        self.panel.box.Add(wx.Button(self.panel, 3, 'Plot'), 0, wx.ALIGN_RIGHT)
        #self.Bind(wx.EVT_BUTTON, self.OnPlot, id=3)
        self.Bind(wx.EVT_BUTTON, lambda evt, values=values: self.OnPlot(evt, values), id=3)
        
        
        #sowas dummes !!!
        
        self.grid.SetRowLabelAlignment(wx.ALIGN_LEFT, wx.ALIGN_CENTRE)
        self.grid.SetRowLabelSize(wx.grid.GRID_AUTOSIZE)        
        self.grid.AutoSize()
        self.panel.box.Layout()
        
    #actual correlation plot
    def OnPlot(self, event, values):
        #prepare values for plot
        print '\nplot pressed\n'
        temp = [[], [], []]
        i = 0
        for x in values[0]:
            if values[2].has_key(x):
                temp[0].append(x)
                temp[1].append(values[2].get(x))
                temp[2].append(values[1][i])
            i = i + 1
        temp.append(values[1][0])
        for l in temp:
            print l
       
        i = 0        
        group = ['RA', 'RARS', 'RAEB', 'healthy']       
        for g in group:
            data1 = []
            data2 = []
            i = 0
            for v in temp[1][:20]:
                if str(v) == str(g):
                    data1.append(temp[2][i])
                i = i + 1
            i = 20         
            for v in temp[1][20:]:
                if str(v) == str(g):
                    data2.append(temp[2][i])
                i = i + 1
            col = 'g'
            if str(g) == 'RA':
                col = 'b'
            if str(g) == 'RARS':
                col = 'pink'
            if str(g) == 'RAEB':
                col = 'r'
            if str(g) == 'healthy':
                col = 'g'
            plt.plot(data1, data2, 'ro', c=col, label=g)
            
            print '\n'
            print data1
            print data2
            
        #plt.legend(bbox_to_anchor=(0., 1.02, 1., .102),loc=3,ncol=2, mode="expand", borderaxespad=0.)
        
        #plt.savefig('samplefigure', bbox_extra_artists=(bbox_extra_artists=(,), bbox_inches='tight',), bbox_inches='tight')
        plt.xlabel('Methylierung')
        plt.ylabel('Genexpression')
        plt.title(values[1][0])
        plt.legend()
        plt.show()                
         
        
        

class DemoPanel(wx.Panel):
    def __init__(self, parent, *args, **kwargs):
        wx.Panel.__init__(self, parent, *args, **kwargs)
        self.parent = parent
        self.box = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.box)
        self.genefield=None
        self.SetAutoLayout(1)
        self.Layout()
        
#not in use
class ChooseMessageDialog(wx.MessageDialog):
    def __init__(self, parent, *args, **kwargs):
        md = wx.MessageDialog.__init__(self, parent, *args, **kwargs)
        pane = wx.Panel(self)
        
        
class MainWindow(wx.Frame):
    def __init__(self, parent, title):
        wx.Frame.__init__(self, parent, title=title, size=(600, 450))           
        self.panel = DemoPanel(self) 
        #self.quote = wx.StaticText(self.panel, label="Your quote: ", pos=(-1,40))
        #Menu auslagern als extra Klasse
        filemenu = wx.Menu()
        
        menuOpenFile = filemenu.Append(wx.ID_OPEN, "&Open file", " Open a File")
        filemenu.AppendSeparator()  
              
        menuOpenFilepath = filemenu.Append(wx.ID_PREVIEW_FIRST, "&Open path", " Open a filepath and search for a specific gene")
        filemenu.AppendSeparator()
        
        menuClearFile = filemenu.Append(wx.ID_CLEAR, "&Clear all", " Clear the Text")
        filemenu.AppendSeparator()
        
        menuAbout = filemenu.Append(wx.ID_ABOUT, "&About", " Information about this program")              
        filemenu.AppendSeparator()
        
        menuExit = filemenu.Append(wx.ID_EXIT, "&Exit", " Terminate the program") 
        self.statusbar = self.CreateStatusBar()   
               
        # Events
        self.Bind(wx.EVT_MENU, self.OnAbout, menuAbout)
        self.Bind(wx.EVT_MENU, self.OnOpenFile, menuOpenFile)  
        self.Bind(wx.EVT_MENU, self.OnExit, menuExit)  
        self.Bind(wx.EVT_MENU, self.OnClear, menuClearFile)
        self.Bind(wx.EVT_MENU, self.OnPreview_First,menuOpenFilepath)
        
        # Menubar.
        menuBar = wx.MenuBar()
        menuBar.Append(filemenu, "&File") # Adding the "filemenu" to the MenuBar
        
        self.SetMenuBar(menuBar)  # Adding the MenuBar to the Frame content.
        self.Show(True)
       
    
    def OnOpenFile(self, event):
        """ Open a file"""
        self.dirname = ''
        dlg = wx.FileDialog(self, "Choose a file",self.dirname, "", "*.*", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()

            self.path = os.path.join(self.dirname, self.filename)
            print 'new path'
            print os.path.abspath(self.path)
            print os.path.basename(self.path)
            print os.path.dirname(self.path)
            #self.path = absolutePathtoRelative(self.path)
            #print self.path
            
            if self.path.find('.tsv') != -1:                                    
                dim = countRowsNCols(self.path);
                self.grid = wx.grid.Grid(self.panel)                
                self.grid.CreateGrid(dim[1], dim[0])
                self.panel.box.Add(self.grid, 1, wx.EXPAND)
                self.grid.Bind(wx.grid.EVT_GRID_CELL_RIGHT_CLICK, self.showPopupMenu)
                f = open(self.path, 'r')
                
                head = f.readline().split('\t')
                i = 0
                for v in head:
                    self.grid.SetColLabelValue(i, v)
                    i = i + 1              
                
                l = 0
                for line in f:                    
                    rows = line.split('\t')
                    c = 0
                    for col in rows:
                        self.grid.SetCellValue(l, c, col)
                        c += 1                        
                    l += 1                        
                f.close()
                
                self.grid.SetColMinimalAcceptableWidth(0)   
                self.grid.SetRowLabelSize(wx.grid.GRID_AUTOSIZE) 
                #disable not needed cols   
                
                try:
                    f = open('grid/gridtitles.tsv', 'r')
                    titlehead = f.readline().strip().split('\t')
                    f.close()

                    i=0
                    head = '\t'.join(head).strip().split('\t')
                    for v in head:
                        if v not in titlehead:
                            self.grid.SetColSize(i, 0)
                        i=i+1

                    """
                    for i in range(dim[0]):
                        if i not in titles:
                            self.grid.SetColSize(i, 0)
                    """
                    
                except IOError:
                    print 'ERROR ----- Please create \grid\gridtitles.tsv and enter the title names separated by linebreak'
               
                
                self.panel.box.Layout()
                
            else:
                f = open(self.path, 'r')
                self.statusbar.SetLabel(self.path)
                self.quote.SetLabel(f.read())           
                f.close()
        dlg.Destroy()
    
    def OnPreview_First(self,event):
        self.dirname = ''
        dlg = wx.FileDialog(self, "Choose a file", self.dirname, "", "*.*", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            self.path = os.path.join(self.dirname, self.filename)
            #self.path = absolutePathtoRelative(self.path)
            
            if self.path.find('.tsv') != -1:      
                   
                self.statusbar.PushStatusText(self.path)   
                """
                self.hsizer = wx.BoxSizer(wx.HORIZONTAL)
                self.hsizer.Add(wx.Button(self.panel, 7, 'search'),1)
                self.Bind(wx.EVT_BUTTON, self.OnSearch, id=7)
                self.genefield = wx.TextCtrl(self.panel)
                self.hsizer.Add(genefield,3,wx.EXPAND)
                self.panel.box.Add(self.hsizer,0, wx.EXPAND)    
                #self.panel.box.Add(wx.Button(self.panel, 3, 'Plot2')) 
                self.SetSizer(self.panel.box)
                """
                
                self.genefield = wx.TextCtrl(self.panel)
                self.panel.box.Add(self.genefield,0,wx.EXPAND)
                
                self.hsizer = wx.BoxSizer(wx.HORIZONTAL)
                self.hsizer.Add(wx.Button(self.panel, 7, 'search'),1)
                self.hsizer.Add(wx.Button(self.panel, 8, 'reset grid'),1)
                self.panel.box.Add(self.hsizer)
                
                #self.panel.box.Add(wx.Button(self.panel, 7, 'search'),0)
                self.Bind(wx.EVT_BUTTON, self.OnSearch, id=7)
                self.Bind(wx.EVT_BUTTON,self.OnReset,id=8)
                #self.panel.box.Add(wx.Button(self.panel, 3, 'Plot3'),1)
                
                #self.hsizer.Layout()
                self.panel.box.Layout()
                
    def OnReset(self, event):
        #remove the grid from sizer
        self.panel.box.Remove(self.grid)
        for child in self.panel.GetChildren(): 
            if type(child)==wx.grid.Grid:
                child.Destroy() 
        self.panel.box.Layout()
              
    def OnSearch(self,event):
        print self.genefield.GetValue()
        f = open(self.path,'r') 
        head = f.readline().split('\t')      
        i = 0
        for x in head:
            if x.lower().find('gene') != -1:
                genecol = i
                break
            i = i + 1
        f.close()
        print 'OnSerach path '+str(self.path)
        if checkIndexFile(self.path):
            
            f = open(self.path,'r')
            temp = long(getIndexOffset(self.genefield.GetValue(), self.path))
            print 'seek '+str(temp)
            f.seek(temp)
            
            values=[]
            
            
            i = 0
            match = False
            for line in f:
                line = line.split('\t')
                if line[genecol] == self.genefield.GetValue():
                    match = True
                    values.append(line)
                    
                if values and line[genecol] != self.genefield.GetValue() and i>500:
                    print 'Stopped search Gene was found'
                    f.close()
                    break
                if match:
                    i=i+1  
        else:
            values=[]
            f = open(self.path,'r')
            f.readline()
            i = 0
            match = False
            for line in f:
                line = line.split('\t')
                if line[genecol] == self.genefield.GetValue():
                    match = True
                    values.append(line)
                    
                if values and line[genecol] != self.genefield.GetValue() and i>500:
                    print 'Stopped search Gene was found'
                    break
                if match:
                    i=i+1              
            f.close()

        self.grid = wx.grid.Grid(self.panel)                
        self.grid.CreateGrid(len(values),len(head))
        self.panel.box.Add(self.grid, 1, wx.EXPAND)
        self.grid.Bind(wx.grid.EVT_GRID_CELL_RIGHT_CLICK, self.showPopupMenu)
        
        i = 0
        for v in head:
            self.grid.SetColLabelValue(i, v)
            i = i + 1              
                
        for i in range(len(values)):
            for j in range(len(head)):          
                self.grid.SetCellValue(i, j,values[i][j])
   
                
        self.grid.SetColMinimalAcceptableWidth(0)   
        self.grid.SetRowLabelSize(wx.grid.GRID_AUTOSIZE) 
        self.panel.box.Layout()     
    
    def OnAbout(self, event):
        dlg = wx.MessageDialog(self, "A small text editor", "About Sample Editor", wx.OK)
        dlg.ShowModal() # Show it
        dlg.Destroy() # finally destroy it when finished.
    
    def OnClear(self, event):
        #remove the grid from sizer
        self.panel.box.Remove(self.grid)
        self.panel.box.Layout()
        #remove all childs from the panel !
        for child in self.panel.GetChildren(): 
            child.Destroy() 
        self.panel.box.Layout()
        
    
    def OnExit(self, event):
        print 'exit'
        self.Close()
        
    #correlation plot    
    def OnPreview(self, event):
        #GetNumberCols
        row = self.grid.GetSelectedRows()[0]
        cols = self.grid.GetNumberCols()
        values = [[]]
        values.append([])
        #Hier werden die ausgewaehlten Werte und die Spaltennamen fuer die Weitergabe zusammengefasst
        for x in range(cols):
            values[0].append(str(self.grid.GetColLabelValue(x)))
            values[1].append(str(self.grid.GetCellValue(row, x)))

        #ChooseFrame(self, values=values).Show()

        mdict = getGroupMapping("maps/KorrelationMapping.csv")

        values.append(mdict)
        for v in values:
            print v

        pos = [[],[]]

        for k in mdict:
            pos[0].append(k)
            temp = []
            for v in mdict.get(k):
                temp.append(values[0].index(v))
            pos[1].append(temp)

        print pos

        #g[0] = ['healthy', 'pre-5-Aza', 'post-5-Aza', 'post-Len', 'pre-Len']
        pos[0][4]='MDS low+int-1'
        pos[0][1]='MDS int-2+high'

        fig = plt.figure()
        ax = plt.subplot(111)

        i = 0
        for g in pos[0]:
            meth = []
            exp = []
            for v in pos[1][i]:
                if v < 78:
                    meth.append(float(values[1][v]))
                elif v > 78:
                    exp.append(float(values[1][v]))
                else:
                    print 'MappingERROR'
            i = i+1

            col = 'g'
            if str(g) == 'healthy':
                col = 'g'
            if str(g) == 'MDS int-2+high':
                col = 'r'
            if str(g) == 'post-5-Aza':
                col = 'y'
            if str(g) == 'post-Len':
                col = 'c'
            if str(g) == 'MDS low+int-1':
                col = 'b'

            ax.plot(meth, exp, 'ro', c=col, label=g)


            print '\n'
            print g
            print meth
            print exp

        #plt.legend(bbox_to_anchor=(0., 1.02, 1., .102),loc=3,ncol=2, mode="expand", borderaxespad=0.)


        #plt.savefig('samplefigure', bbox_extra_artists=(bbox_extra_artists=(,), bbox_inches='tight',), bbox_inches='tight')


        plt.xlabel('Methylierung')
        plt.ylabel('Genexpression')
        plt.title(values[1][0])

        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])

        # Put a legend to the right of the current axis
        #ax.legend(loc='upper left', bbox_to_anchor=(1, 0.5))
        ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.show()





    #splicing plot   
    def OnPrint(self, event):
        #change this if first col is gene names !!!
        row = self.grid.GetSelectedRows()[0]
        cols = self.grid.GetNumberCols()
        values = [[]]
        #values.append([])
        #Hier werden die ausgewaehlten Werte und die Spaltennamen fuer die Weitergabe zusammengefasst
        for x in range(cols):
            values[0].append(str(self.grid.GetColLabelValue(x)))
        #find col index with gene in title
        i = 0
        for x in values[0]:
            if x.lower().find('gene') != -1:
                genecol = i
                break
            i = i + 1
            
        #zeilen suchen in der das gewaehlte gen vorkommt
        genename = str(self.grid.GetCellValue(row, genecol))
        
        if row == 0:
            while str(self.grid.GetCellValue(row, genecol)) == genename:        
                temp = []
                for x in range(cols):
                    temp.append(str(self.grid.GetCellValue(row, x)))     
                values.append(temp)
                row = row + 1        
        else:     
            while str(self.grid.GetCellValue(row, genecol)) == genename and str(self.grid.GetCellValue(row, genecol)) != None :
                #die null abfragen waren weg
                row = row - 1
                if row<0:
                    break
            row = row + 1
            try:
                while str(self.grid.GetCellValue(row, genecol)) == genename:     
                    temp = []
                    for x in range(cols):
                        temp.append(str(self.grid.GetCellValue(row, x)))     
                    values.append(temp)
                    row = row + 1
            except AssertionError:
                print 'error'                
                #if row >= self.grid.GetNumberRows():break
                  
        mdict = getGroupMapping("maps/groupmapping.csv")
        values.append(mdict)

        for v in values:
            print v
            
        grouppos = []
        for g in range(len(set(mdict.values()))):
            grouppos.append([])
        i = 0   
        for x in values[0]:
            if values[len(values) - 1].has_key(x):
                #if values[len(values)-1].get(x)=='RA':
                if(values[len(values) - 1].get(x)) == 'RA':
                    grouppos[0].append(i)
                elif(values[len(values) - 1].get(x)) == 'RARS':
                    grouppos[1].append(i)
                elif(values[len(values) - 1].get(x)) == 'RAEB':
                    grouppos[2].append(i)
                elif(values[len(values) - 1].get(x)) ==  'healthy':
                    grouppos[3].append(i)
                else:
                    print 'Fehler bei Gruppenpositionfindung'
                      
            i = i + 1
        
        drawarray=[]
        for g in range(len(set(mdict.values()))):
            drawarray.append([])
        for line in values[1:len(values)-1]:
            ra = []
            rars=[]
            raeb=[]
            healthy=[]
            for g in grouppos[0]:
                ra.append(line[g])
            drawarray[0].append(ra)
            
            for g in grouppos[1]:
                rars.append(line[g])
            drawarray[1].append(rars)
            
            for g in grouppos[2]:
                raeb.append(line[g])
            drawarray[2].append(raeb)
            
            for g in grouppos[3]:
                healthy.append(line[g])
            drawarray[3].append(healthy)
        #Mittelwerte und Standardabweichung berechnen
        meandrawarray = deepcopy(drawarray)
        devdrawarray = deepcopy(drawarray)
        for i in range(len(drawarray)):
            for j in range(len(drawarray[i])):
                m = map(float, meandrawarray[i][j])
                drawarray[i][j]=m
                meandrawarray[i][j]=mean(m)
                devdrawarray[i][j]=std(m)
                
        
        print '\nrohwerte'
        for x in drawarray:
            print x
        print '\nmittelwerte'        
        for v in meandrawarray:
            print v
        print '\nstandardabweichung'
        for v in devdrawarray:
            print v
        print '\n'
                       
        exonintron = getExon(genecol, values[1:len(values) - 1])
        exon=exonintron[0]
        intron=exonintron[1]        
     
        fig = plt.figure(1)
        #plt.subplot(211)
        
        gs = gridspec.GridSpec(2,1,height_ratios=[1,10*(math.exp( -0.3*len(exon) ))+2])
        
        ax3 = fig.add_subplot(gs[0])


        print '\nexon'
        for e in exon:print e
        print 'intron'
        for i in intron:print i
        
        j=len(exon)        
        for a in exon:
            i=0
            done = False
            while i < len(a)-1:
                if not done:
                    if a[i][0]+1 == a[i+1][0] and a[i][1] == a[i+1][1]:
                        s=i
                        while a[i][0]+1 == a[i+1][0] and a[i][1] == a[i+1][1]:
                            i=i+1
                            if i == len(a)-1:
                                #print str(s)+'-'+str(i)
                                rect = patches.Rectangle((a[s][0]-0.5,j), i-s+1, 0.5, edgecolor='black',facecolor='grey')
                                ax3.add_patch(rect)
                                done = True
                                print 'done'
                                break
                        else:
                            #print str(s)+'-'+str(i)  
                            rect = patches.Rectangle((a[s][0]-0.5,j), i-s+1, 0.5, edgecolor='black',facecolor='grey')
                            ax3.add_patch(rect)          
                            i=i+1
                            print 'continue'
                            if i != len(a)-1:
                                continue
                if not done:     
                    #print i 
                    rect= patches.Rectangle((a[i][0]-0.5,j), 1, 0.5, edgecolor='black',facecolor='grey')
                    ax3.add_patch(rect)                       
                    i=i+1
                    if i == len(a)-1:
                        rect= patches.Rectangle((a[i][0]-0.5,j), 1, 0.5, edgecolor='black',facecolor='grey')
                        ax3.add_patch(rect)
            j=j-1
        
        j=len(intron)        
        for a in intron:
            i=0
            done = False
            while i < len(a)-1:
                if not done:
                    if a[i][0]+1 == a[i+1][0] and a[i][1] == a[i+1][1]:
                        s=i
                        while a[i][0]+1 == a[i+1][0] and a[i][1] == a[i+1][1]:
                            i=i+1
                            if i == len(a)-1:
                                #print str(s)+'-'+str(i)
                                rect = patches.Rectangle((a[s][0]-0.5,j), i-s+1, 0.5, edgecolor='black',facecolor='white')
                                ax3.add_patch(rect)
                                done = True
                                print 'done'
                                break
                        else:
                            #print str(s)+'-'+str(i)  
                            rect = patches.Rectangle((a[s][0]-0.5,j), i-s+1, 0.5, edgecolor='black',facecolor='white')
                            ax3.add_patch(rect)          
                            i=i+1
                            print 'continue'
                            if i != len(a)-1:
                                continue
                if not done:     
                    #print i 
                    rect= patches.Rectangle((a[i][0]-0.5,j), 1, 0.5, edgecolor='black',facecolor='white')
                    ax3.add_patch(rect)                       
                    i=i+1
                    if i == len(a)-1:
                        rect= patches.Rectangle((a[i][0]-0.5,j), 1, 0.5, edgecolor='black',facecolor='white')
                        ax3.add_patch(rect)
            j=j-1
            
        """   
            for e in t:
                rect = matplotlib.patches.Rectangle((e[0],i), 1, 0.5, edgecolor='black',facecolor='grey')
                ax.add_patch(rect)
            i=i-1
        """
        
        #plt.xlim([-1,exon[0][len(exon)-1][0]+20])
        #plt.ylim([0,2*len(exon)+5])
        plt.xlim([-1,len(values)-2])
        plt.ylim([0,len(exon)+1])
        plt.ylabel('Transkript')
        plt.grid(True)
        plt.title(values[1][1])
        #x und y achsen steps als ganzzahlen
        plt.xticks(range(len(meandrawarray[0])))
        plt.yticks(range(len(exon)+1))
        
        #ax2 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        fig.add_subplot(gs[1])
        #ax2.grid(b=True, which='major')
        rar = range(len(meandrawarray[0]))
        rarsr = range(len(meandrawarray[1]))
        raebr = range(len(meandrawarray[2]))
        healthyr = range(len(meandrawarray[3]))
        ra = meandrawarray[0]
        rars = meandrawarray[1]
        raeb = meandrawarray[2]
        healthy = meandrawarray[3]
        
        rac,rarsc,raebc,healthyc = 'b','g','r','c'
        
        #mittelwert durch kurven plotten
        plt.plot(rar,ra,color=rac,linewidth=2,label='RA')
        plt.plot(rarsr,rars,color=rarsc,linewidth=2,label='RARS')
        plt.plot(raebr,raeb,color=raebc,linewidth=2,label='RAEB')
        plt.plot(healthyr,healthy,color=healthyc,linewidth=2,label='healthy')
        plt.legend()
        
        
        i=0
        while i<len(rar):
            rar[i]=rar[i]-0.1
            rarsr[i]=rarsr[i]-0.05
            raebr[i]=raebr[i]+0.05
            healthyr[i]=healthyr[i]+0.1
            i=i+1  
        
        #Errorbar durch standardabweichung zeichnen fmt blockiert die farb belegung
        plt.errorbar(range(len(ra)), ra, yerr=devdrawarray[0],fmt=None,color=rac,capthick=2,label="series 1",capsize=5)
        #print 'laenge w : '+ str(len(w)) +' w :'+str(w)
        plt.errorbar(range(len(rars)), rars, yerr=devdrawarray[1],fmt=None,color=rarsc,capthick=2,label="series 1",capsize=5)
        plt.errorbar(range(len(raeb)), raeb, yerr=devdrawarray[2],fmt=None,color=raebc,capthick=2,label="series 1",capsize=5)
        plt.errorbar(range(len(healthy)), healthy, yerr=devdrawarray[3],fmt=None,color=healthyc,capthick=2,label="series 1",capsize=5)
        
        """ 
        i=0
        for x in drawarray:
            i=0
            for y in x:
                for z in y:
                    plt.plot(i,z,'b.')
                i=i+1
        """
        #punkte ploten
        g=0
        for x in drawarray:
            try :
                i=0
                while True:
                    temp=[]
                    for y in x:
                        temp.append(y[i])
                    #print 'from group '+str(g)+' pop '+str(i)+' '+ str(temp)
                    
                    col='b.'
                    if g==0:
                        col=rac+'.'
                        dimr=rar
                    elif g==1:
                        col=rarsc+'.'
                        dimr=rarsr
                    elif g==2:
                        col=raebc+'.'
                        dimr=raebr
                    elif g==3:
                        col=healthyc+'.'
                        dimr=healthyr
                    #print 'laenge temp : '+str(len(temp))+' temp '+str(temp)

                    plt.plot(dimr,temp,col)
                
                    i=i+1
            
            except IndexError:
                pass
            g=g+1
      
        maxi = 1
        for x in meandrawarray:
                maxi = max(x)
        mini = 0
        for x in meandrawarray:
                mini=min(x)
                
        plt.ylim((mini-(mini/5)), maxi+(maxi/5))   
        plt.xlim(-1,len(values)-2)   
        plt.ylabel('Genexpression in AU')
        plt.xlabel('Probe Set')
        plt.grid(True)
        plt.xticks(range(len(meandrawarray[0])))

        plt.show()       
        
        
    def showPopupMenu(self, event):
        self.grid.SelectRow(event.GetRow())
        print event.GetRow()
        
        """ TODO statusbar sollte auch beim popupmenu mouseover zeigen"""      
        if not hasattr(self, "popupID1"):
            self.popupID1 = wx.NewId()
            self.popupID2 = wx.NewId()
            self.popupID3 = wx.NewId()
            # make a menu
 
        popmenu = wx.Menu()
        # Show how to put an icon in the menu
        menuPlot2 = popmenu.Append(wx.ID_PREVIEW, "&Korrelationsplot", "Draw Correlation Plot")
        menuPlot3 = popmenu.Append(wx.ID_PRINT, "&Genexpressions- und Transkriptplot", "Draw Splicing Plot")
        #popmenu.Append(self.popupID3, "Splicing Plot")
        
        
        #bindings for popmenu
        self.Bind(wx.EVT_MENU, self.OnPreview, menuPlot2) 
        self.Bind(wx.EVT_MENU, self.OnPrint, menuPlot3) 
        # Popup the menu.  If an item is selected then its handler
        # will be called before PopupMenu returns.
        self.PopupMenu(popmenu)
        popmenu.Destroy()
        

app = wx.App(False)
frame = MainWindow(None, "")
app.MainLoop()
