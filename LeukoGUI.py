'''
Created on 23.04.2013
@author: Axel Wilbertz
'''

import wx
import os
import wx.grid
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from copy import deepcopy
from numpy import mean,std
import File
import Plots
import math 

class DemoPanel(wx.Panel):
    def __init__(self, parent, *args, **kwargs):
        wx.Panel.__init__(self, parent, *args, **kwargs)
        self.parent = parent
        self.box = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.box)
        self.genefield=None
        self.SetAutoLayout(1)
        self.box.Fit(self)
        self.Layout()

     
class MainWindow(wx.Frame):
    
    fileobj = None
    grid = None
    
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
        self.Bind(wx.EVT_MENU, self.OnClearAll, menuClearFile)
        self.Bind(wx.EVT_MENU, self.OnOpenFilepath,menuOpenFilepath)
        
        # Menubar.
        
        opmenu = wx.Menu()
        menuCalcLogLinear = opmenu.Append(wx.ID_ANY, "&Log/Linear Calculation", "Do some awesome calculation")
        opmenu.AppendSeparator()  
        
        
        self.Bind(wx.EVT_MENU, self.OnCalcLogLinear,menuCalcLogLinear)
        
        
        menuBar = wx.MenuBar()
        menuBar.Append(filemenu, "&File") # Adding the "filemenu" to the MenuBar
        menuBar.Append(opmenu, "&Fileoperations") 
        
        self.SetMenuBar(menuBar)  # Adding the MenuBar to the Frame content.
        self.Show(True)
    
    def OnCalcLogLinear(self,event):   
        hsizer = wx.BoxSizer(wx.HORIZONTAL)
        browserButton = wx.Button(self.panel,22, 'browse')
        adressField = wx.TextCtrl(self.panel)
        
        hsizer.Add(browserButton,1,wx.EXPAND)
        hsizer.Add(adressField,1,wx.EXPAND)
        newAdressField = wx.TextCtrl(self.panel)
        hsizer.Add(wx.StaticText(self.panel,-1,'New File: ',style=wx.ALIGN_CENTRE),1,wx.EXPAND)
        hsizer.Add(newAdressField,1,wx.EXPAND)
        
        
        #self.panel.box.Add(hsizer2)            
        self.panel.box.Add(hsizer,wx.EXPAND)
        
        #hsizer2 = wx.BoxSizer(wx.HORIZONTAL)
        
     
        self.panel.box.Layout()
        
        print 'Log Lin'
    
    def OnOpenFile(self, event):
        """ Open a file"""
        self.dirname = ''
        dlg = wx.FileDialog(self, "Choose a file",self.dirname, "", "*.*", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()

            #self.path = os.path.join(self.dirname, self.filename)
            self.fileobj = File.FileObject(os.path.join(self.dirname, self.filename),'\t')
            self.fileobj.getTitles()   
            col = self.fileobj.findGeneColumn("maps/GeneColumn.tsv")
            self.fileobj.getGenePos(col)
       
            if self.fileobj.filepath.find('.tsv') != -1:           
                                         
                dim = self.fileobj.countRowsNCols(self.fileobj.filepath)

                self.grid = wx.grid.Grid(self.panel)                
                self.grid.CreateGrid(dim[1], dim[0])
                self.panel.box.Add(self.grid, 1, wx.EXPAND)
                self.grid.Bind(wx.grid.EVT_GRID_CELL_RIGHT_CLICK, self.showPopupMenu)
                
                             
                
                i = 0
                for v in self.fileobj.titles:
                    self.grid.SetColLabelValue(i, v)
                    i = i + 1      
                            
                self.fileobj.getAllValues()                    
                
                l = 0
                for row in self.fileobj.values:
                    c = 0
                    for col in row:
                        self.grid.SetCellValue(l, c, col)
                        c += 1                        
                    l += 1     
                
                self.grid.SetColMinimalAcceptableWidth(0)   
                self.grid.SetRowLabelSize(wx.grid.GRID_AUTOSIZE) 
                
                
                #show only needed cols   
                """
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

                except IOError:
                    print 'ERROR ----- Please create \grid\gridtitles.tsv and enter the title names separated by linebreak'
               """
                
                self.panel.box.Layout()
                
            else:
                raise IOError('File has wrong ending. Wrong data structure could be possible. File should be .tsv and tab stop delimitered')
        dlg.Destroy()
        
    #open filepath
    def OnOpenFilepath(self,event):
        self.dirname = ''
        dlg = wx.FileDialog(self, "Choose a file", self.dirname, "", "*.*", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
              
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            #self.path = os.path.join(self.dirname, self.filename)
            #self.path = absolutePathtoRelative(self.path)
            self.fileobj = File.FileObject(os.path.join(self.dirname, self.filename),'\t')          
            self.fileobj.getTitles()   
            col = self.fileobj.findGeneColumn("maps/GeneColumn.tsv")
            self.fileobj.getGenePos(col)
            
            
            if self.fileobj.filepath.find('.tsv') != -1:      
                   
                self.statusbar.PushStatusText(self.fileobj.filepath)   
                
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
            else:
                raise IOError('File has wrong ending. Wrong data structure could be possible. File should be .tsv and tab stop delimitered')
                        
    def OnReset(self, event):
        """
        Reset Searchfield and results displayed in the grid
        """
        self.panel.box.Remove(self.grid)
        for child in self.panel.GetChildren(): 
            if type(child)==wx.grid.Grid:
                child.Destroy() 
        self.panel.box.Layout()
        #del self.fileobj
              
    def OnSearch(self,event):
        print self.genefield.GetValue()
        #genecol
        if self.grid:
            self.OnReset(event)
        
        self.fileobj.getTitles()
        
        fileidx = File.IndexFile(self.fileobj.filepath,self.fileobj.genecolumn)
        
        print 'OnSepostAzach path '+str(self.fileobj.filepath)
        
        if fileidx.checkIndexFile():            
            #searches for the gene in the indexfile
            offset = long(fileidx.getIndexOffset(self.genefield.GetValue()))
            self.fileobj.getOffsetValues(offset,self.genefield.GetValue())
          
        else:
            """
            values=[]
            f = open(self.path,'r')
            f.readline()
            i = 0
            match = False
            for line in f:
                line = line.split('\t')
                if line[self.fileobj.genecolumn] == self.genefield.GetValue():
                    match = True
                    values.append(line)
                    
                if values and line[self.fileobj.genecolumn] != self.genefield.GetValue() and i>500:
                    print 'Stopped search Gene was found'
                    break
                if match:
                    i=i+1              
            f.close()
            """
            raise ValueError('CheckIndexFile False - Never the case because indexfile creates itself inside boolean call')        
        
        self.grid = wx.grid.Grid(self.panel)                
        self.grid.CreateGrid(len(self.fileobj.values),len(self.fileobj.titles))
        self.panel.box.Add(self.grid, 1, wx.EXPAND)
        self.grid.Bind(wx.grid.EVT_GRID_CELL_RIGHT_CLICK, self.showPopupMenu)
        
        i = 0        
        for t in self.fileobj.titles:
            self.grid.SetColLabelValue(i, t)
            i = i + 1              
                
        for i in range(len(self.fileobj.values)):
            for j in range(len(self.fileobj.titles)):          
                self.grid.SetCellValue(i, j,self.fileobj.values[i][j])   
                
        self.grid.SetColMinimalAcceptableWidth(0)   
        self.grid.SetRowLabelSize(wx.grid.GRID_AUTOSIZE) 
        self.panel.box.Layout()     
    
    def OnAbout(self, event):
        dlg = wx.MessageDialog(self, "A small text editor", "About Sample Editor", wx.OK)
        dlg.ShowModal() # Show it
        dlg.Destroy() # finally destroy it when finished.
    
    def OnClearAll(self, event):
        """
        reset gui and loaded filepaths. reset file object
        """
        if self.grid:
            self.panel.box.Remove(self.grid)
        #remove all childs from the panel !
        for child in self.panel.GetChildren(): 
            child.Destroy() 
        self.panel.box.Layout()
        del self.fileobj
        

    def OnExit(self, event):
        print 'exit'
        self.Close()
        
    #correlation plot    
    def OnCorrelationPlot(self, event):
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
       
        titles = values[0]
        values = values[1]
        genecolumn = 0
        
        
        mapper = File.Mapper('maps\KorrelationMapping.csv',':');
        mapper.extractGroupMapping()
        print 'titles: '+str(titles)
        print 'values: '+str(values)
        print 'mapper: '+str(mapper.mapperdict)
        mapper.findPositions(titles,mapper.mapperdict)
        print 'pos: '+str(mapper.positions)
        corplot = Plots.CorrelationPlot(titles,values,genecolumn,['g','r','y','c','b'],mapper.positions);
        corplot.plotValues()

    #splicing plot   
    def OnSplicingPlot(self, event):
        #change this if first col is gene names !!!
        row = self.grid.GetSelectedRows()[0]
        cols = self.grid.GetNumberCols()
        values = [[]]
        #values.append([])
        #Hier werden die ausgewaehlten Werte und die Spaltennamen fuer die Weitergabe zusammengefasst
        for x in range(cols):
            values[0].append(str(self.grid.GetColLabelValue(x).strip()))
        #find col index with gene in title
        
            
        #zeilen suchen in der das gewaehlte gen vorkommt
        genename = str(self.grid.GetCellValue(row, self.fileobj.genecolumn))
        
        if row == 0:
            while str(self.grid.GetCellValue(row, self.fileobj.genecolumn)) == genename:        
                temp = []
                for x in range(cols):
                    temp.append(str(self.grid.GetCellValue(row, x)))     
                values.append(temp)
                row = row + 1        
        else:     
            while str(self.grid.GetCellValue(row, self.fileobj.genecolumn)) == genename and str(self.grid.GetCellValue(row, self.fileobj.genecolumn)) != None :
                #die null abfragen waren weg
                row = row - 1
                if row<0:
                    break
            row = row + 1
            try:
                while str(self.grid.GetCellValue(row, self.fileobj.genecolumn)) == genename:     
                    temp = []
                    for x in range(cols):
                        temp.append(str(self.grid.GetCellValue(row, x)))     
                    values.append(temp)
                    row = row + 1
            except AssertionError:
                print 'Assertion error'                
                #if row >= self.grid.GetNumberRows():break
        
        titles = values[0]
        values = values[1:]
        genecolumn = 1
                
        mapper2 = File.Mapper('maps\SplicingTranscriptExpressionMapping.csv',':');
        mapper2.extractGroupMapping()
        mapper2.findPositions(titles,mapper2.mapperdict)
        
        file3 = File.RefGene('refgene/refGene.tsv','indexfiles/refGene_index.tsv',genecolumn,'\t',12,9,10,7,8,titles,values)
        file3.extractRefGeneIndexFileEntries()
        file3.extractRefGeneValues()
        file3.calcIntronExon()    
               
        splicingplot= Plots.SplicingPlot(titles,values,self.fileobj.genecolumn,['b','g','r','c','k'],mapper2.positions,file3.intron,file3.exon)
        splicingplot.setMDSSubtypes()
        splicingplot.calcMeanANDStd()
        splicingplot.plotValues()
         
        
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
        menuCorItem = popmenu.Append(wx.ID_PREVIEW, "&Korrelationsplot", "Draw Correlation Plot")
        menuSplicingItem = popmenu.Append(wx.ID_PRINT, "&Genexpressions- und Transkriptplot", "Draw Splicing Plot")
        #popmenu.Append(self.popupID3, "Splicing Plot")
        
        
        #bindings for popmenu
        self.Bind(wx.EVT_MENU, self.OnCorrelationPlot, menuCorItem) 
        self.Bind(wx.EVT_MENU, self.OnSplicingPlot, menuSplicingItem) 
        # Popup the menu.  If an item is selected then its handler
        # will be called before PopupMenu returns.
        self.PopupMenu(popmenu)
        popmenu.Destroy()
        

app = wx.App(False)
frame = MainWindow(None, "")
app.MainLoop()
