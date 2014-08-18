'''
Created on 18.07.2014

@author: Axel Wilbertz
'''
import matplotlib.pyplot as plt
from copy import deepcopy
from numpy import mean,std
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
import math

class Plots(object):
    
    titles=[];
    values=[];
    colours=[];
    pos=[[],[]];
    genecolumn = -1;

    def __init__(self,titles,values,genecolumn,colours,pos):
        '''
        Constructor
        '''
        self.genecolumn=genecolumn
        self.titles=titles;
        self.values = values;
        #clours =['g','r','y','c','b']
        self.colours = colours;
        self.pos = pos;
        

class SplicingPlot(Plots):
    """
    Plots the Methylation and Expression of a Gene through Matplotlib
    """
    
    drawarray = []
    mean = []
    std = []
    
    healthy = []
    postLen=[]
    MDShigh=[]
    postAza=[]
    MDSlow=[]
        
    exon= []
    intron = []
    
    def __init__(self,titles,values,genecolumn,colours,pos,intron,exon):
        '''
        Constructor
        '''
        self.genecolumn=genecolumn
        self.titles=titles;
        self.values = values;
        #clours =['g','r','y','c','b']
        self.colours = colours;
        self.pos = pos;
        self.intron=intron
        self.exon=exon
   
    
    def plotValues(self):       
        
        fig = plt.figure(1)
        
        gs = gridspec.GridSpec(2,1,height_ratios=[1,10*(math.exp( -0.3*len(self.exon) ))+2])      
        
        ax3 = fig.add_subplot(gs[0])
    
    
        j=len(self.exon)        
        for a in self.exon:
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
                                rect = patches.Rectangle((a[s][0]-0.25,j), i-s+0.5, 0.5, edgecolor='black',facecolor='grey')
                                ax3.add_patch(rect)
                                done = True
                                break
                        else:
                            #print str(s)+'-'+str(i)  
                            rect = patches.Rectangle((a[s][0]-0.25,j), i-s+0.5, 0.5, edgecolor='black',facecolor='grey')
                            ax3.add_patch(rect)          
                            i=i+1
                            if i != len(a)-1:
                                continue
                if not done:     
                    #print i 
                    rect= patches.Rectangle((a[i][0]-0.25,j), 0.5, 0.5, edgecolor='black',facecolor='grey')
                    ax3.add_patch(rect)                       
                    i=i+1
                    if i == len(a)-1:
                        rect= patches.Rectangle((a[i][0]-0.25,j), 0.5, 0.5, edgecolor='black',facecolor='grey')
                        ax3.add_patch(rect)
            j=j-1
        
        
        
        #plt.xlim([-1,exon[0][len(exon)-1][0]+20])
        #plt.ylim([0,2*len(exon)+5])
        plt.xlim([-1,len(self.values)])
        plt.ylim([0,len(self.exon)+1])
        plt.ylabel('Transkript')
        plt.grid(True)
        plt.title(self.values[1][1])
        #x und y achsen steps als ganzzahlen
        plt.xticks(range(len(self.mean[0])))
        plt.yticks(range(len(self.exon)+1))
        
        box = ax3.get_position()
        ax3.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
        # Put a legend to the right of the current axis
        #ax.legend(loc='upper left', bbox_to_anchor=(1, 0.5))
        
        #ax2 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        ax2 = fig.add_subplot(gs[1])
        #ax2.grid(b=True, which='major')
        
            
        healthyr = range(len(self.mean[0]))
        postLenr = range(len(self.mean[1]))
        MDShighr = range(len(self.mean[2]))
        postAzar = range(len(self.mean[3]))
        MDSlowr = range(len(self.mean[4]))
        
        self.healthy = self.mean[0]
        self.postLen = self.mean[1]
        self.MDShigh = self.mean[2]
        self.postAza = self.mean[3]
        self.MDSlow = self.mean[4]
        
        postAzac,healthyc,postLenc,MDShighc,MDSlowc = self.colours
        
        #mittelwert durch kurven plotten
        p1 = plt.plot(postAzar,self.postAza,color=postAzac,linewidth=2,label='postAza')
        p2 = plt.plot(postLenr,self.postLen,color=postLenc,linewidth=2,label='postLen')
        p3 = plt.plot(MDShighr,self.MDShigh,color=MDShighc,linewidth=2,label='MDShigh')
        p4 = plt.plot(MDSlowr,self.MDSlow,color=MDSlowc,linewidth=2,label='MDSlow')
        p5 = plt.plot(healthyr,self.healthy,color=healthyc,linewidth=2,label='healthy')
        
        
        plt.legend()
        box = ax2.get_position()
        ax2.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
        ax2.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        
        
        i=0
        while i<len(healthyr):
            healthyr[i]=healthyr[i]-0.1
            postLenr[i]=postLenr[i]-0.05
            MDShighr[i]=MDShighr[i]+0.05
            postAzar[i]=postAzar[i]+0.1
            MDSlowr[i]=MDSlowr[i]+0.15            
            i=i+1  
        
        #Errorbar durch standardabweichung zeichnen fmt blockiert die farb belegung
        plt.errorbar(range(len(self.healthy)), self.healthy, yerr=self.std[0],color=healthyc,capthick=2,label="series 1",capsize=5)
        #print 'laenge w : '+ str(len(w)) +' w :'+str(w)
        plt.errorbar(range(len(self.postLen)), self.postLen, yerr=self.std[1],color=postLenc,capthick=2,label="series 1",capsize=5)
        plt.errorbar(range(len(self.MDShigh)), self.MDShigh, yerr=self.std[2],color=MDShighc,capthick=2,label="series 1",capsize=5)
        plt.errorbar(range(len(self.postAza)), self.postAza, yerr=self.std[3],color=postAzac,capthick=2,label="series 1",capsize=5)
        plt.errorbar(range(len(self.MDSlow)), self.MDSlow, yerr=self.std[4],color=MDSlowc,capthick=2,label="series 1",capsize=5)
        
        #punkte ploten
        g=0
        for x in self.drawarray:
            try :
                i=0
                while True:
                    temp=[]
                    for y in x:
                        temp.append(y[i])
                    #print 'from group '+str(g)+' pop '+str(i)+' '+ str(temp)
                    
                    col='b.'
                    if g==0:
                        col=healthyc+'.'
                        dimr=healthyr
                    elif g==1:
                        col=postLenc+'.'
                        dimr=postLenr
                    elif g==2:
                        col=MDShighc+'.'
                        dimr=MDShighr
                    elif g==3:
                        col=postAzac+'.'
                        dimr=postAzar
                    elif g==4:
                        col=MDSlowc+'.'
                        dimr=MDSlowr
                    
                    #print 'laenge temp : '+str(len(temp))+' temp '+str(temp)                  
                    plt.plot(dimr,temp,col)
                
                    i=i+1
            
            except IndexError:
                pass
            g=g+1
      
        maxi = 1
        for x in self.mean:
                maxi = max(x)
        mini = 0
        for x in self.mean:
                mini=min(x)
                
        
        plt.ylim((mini-(mini/5)), maxi+(maxi/5))   
        plt.xlim(-1,len(self.values))   
        plt.ylabel('Genexpression in AU')
        plt.xlabel('Probe Set')
        plt.grid(True)
        plt.xticks(range(len(self.mean[0])))
        
        plt.show()       
        
    def setMDSSubtypes(self):
        self.drawarray=None
        self.drawarray=[]        
                
        for g in range(len(self.pos[0])):
            self.drawarray.append([])
        
        for line in self.values:
            
            healthy = []
            postLen=[]
            MDShigh=[]
            postAza=[]
            MDSlow=[]
            
            """
            change this !!! first row mds subtypes like healthy second row to last row values for subtypes
            """
            
            for g in self.pos[1]:
                healthy.append(line[g])
            self.drawarray[0].append(healthy)
            
            for g in self.pos[2]:
                postLen.append(line[g])
            self.drawarray[1].append(postLen)
            
            for g in self.pos[3]:
                MDShigh.append(line[g])
            self.drawarray[2].append(MDShigh)
            
            for g in self.pos[4]:
                postAza.append(line[g])
            self.drawarray[3].append(postAza)
            
            for g in self.pos[5]:
                MDSlow.append(line[g])
            self.drawarray[4].append(MDSlow)
            
            self.healthy = healthy
            self.postLen = postLen
            self.MDShigh = MDShigh
            self.postAza = postAza
            self.MDSlow = MDSlow            
            
        return [healthy,postLen,MDShigh,postAza,MDSlow]
          
    def calcMeanANDStd(self):
        #Mittelwerte und Standardabweichung berechnen
        meandrawarray = deepcopy(self.drawarray)
        devdrawarray = deepcopy(self.drawarray)
        for i in range(len(self.drawarray)):
            for j in range(len(self.drawarray[i])):
                m = map(float, meandrawarray[i][j])
                self.drawarray[i][j]=m
                meandrawarray[i][j]=round(mean(m),4)
                devdrawarray[i][j]=round(std(m),4)
               
        self.mean = meandrawarray
        self.std = devdrawarray  
        print '\nmean'
        for m in self.mean:print m
        print 'std dev'
        for s in self.std:print s
        return [meandrawarray,devdrawarray]

class CorrelationPlot(Plots):
    """
    Plots the Expression of a Gene through Matplotlib. Calculates Introns and Exons.
    """
    
    def plotValues(self): 
            fig = plt.figure()
            ax = plt.subplot(111)
    
            i = 0
            for g in self.pos[0]:
                meth = []
                exp = []
                for v in self.pos[1][i]:
                    if v < 78:
                        # meth.append(float(self.values[1][v]))
                        meth.append(float(self.values[v]))
                    elif v > 78:
                        # exp.append(float(self.values[1][v]))
                        exp.append(float(self.values[v]))
                    else:
                        print 'MappingERROR'
                
                col = self.colours[i]
                
                i = i+1
                
                ax.plot(meth, exp, 'ro', c=col, label=g)
        
                print '\n'
                print g
                print meth
                print exp
    
            #plt.legend(bbox_to_anchor=(0., 1.02, 1., .102),loc=3,ncol=2, mode="expand", borderaxespad=0.)
    
            #plt.savefig('samplefigure', bbox_extra_artists=(bbox_extra_artists=(,), bbox_inches='tight',), bbox_inches='tight')
       
            plt.xlabel('Methylierung')
            plt.ylabel('Genexpression')
            plt.title(self.values[self.genecolumn])
    
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
    
            # Put a legend to the right of the current axis
            #ax.legend(loc='upper left', bbox_to_anchor=(1, 0.5))
            ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            plt.show()



