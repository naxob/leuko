import sets
import os
import SortFile
import numpy as np
from scipy import stats
import File



'''
Created on 11.12.2013
New Microarray files added
@author: Axel Wilbertz
'''
from test.test_argparse import WFile
"""
@return: returns index of gen column.if not found return -1
"""
def getgenpos(array, gencolname,*delimiter):      
        if type(array) == str and array.index('.'):
            file = open(array,'r')            
            array = file.readline()
            array = array.split(delimiter[0])   
            file.close()         
        if(type(gencolname) == str):
            if type(array[0])==list:
                return array[0].index(gencolname)
            else:
                return array.index(gencolname)
        elif type(gencolname) == list :
            indexarray = []
            if type(array[0])==list:            
                for elem in gencolname:           
                    indexarray.append(array[0].index(elem))
            else:
                for elem in gencolname:           
                    indexarray.append(array.index(elem))
            return indexarray  
            
def findvalinarray(array,val):
    pos = []
    i = 0
    for x in array:
        if x.find(val) >= 0:
            pos.append(i)
        i+=1
    return pos
        
def countlines(file):
    file.seek(0)
    i=0
    for line in file:
        i+=1
    file.seek(0)
    return i

"""
@param linecount:set 0 for end of file 
"""
def delete(array, col, it):
    temparray = []
    for line in array:
        if line[col] != it:
            temparray.append(line)
    return temparray

def deleteinstring(array, col, it):
    print 'deleting in '+str(array[0][col])
    for line in array:       
        line[col] = line[col].translate(None, it)
    return array

def replaceinarray(array,col,it):    
    for x in col:
        array[x] = array[x].replace(it,'')    
    return array
    
           
def stripper(array, col):
    print '\nstripping array'
    for line in array:
        line[col] = line[col].strip()
    return array

def readfile(filepath, delimiter, linecount):
    print '\nreadfile started...'
    file = open(filepath, 'r')
    temparray = []
    if linecount:        
        i = 1        
        for line in file:       
            temparray.append(line.split(delimiter))
            if i >= linecount:break
            i += 1
    else:
        for line in file:
            temparray.append(line.split(delimiter))
    print '\t' + str(len(temparray)) + ' lines added to array and split by ' + delimiter
    return temparray

def addID(list):  
    i = 0
    for line in list:
        if i == 0:
            line.insert(0, 'ID')
            
        else:
            line.insert(0, str(i - 1))           
        i += 1
    
    return list

def printfile(filepath,linecount,delimiter=0,colpos=0):
    print 'printing '+str(filepath)
    file = open(filepath, 'r')
    if colpos and delimiter:
        i=0
        for line in range(linecount):
            print str(i)+' '+file.readline().split(delimiter)[colpos]
            i+=1
    elif delimiter:
        for line in range(linecount):
            print file.readline().split(delimiter)
    else:
        for line in range(linecount):
            print file.readline().strip('\n')
    
    file.close()
    return

def printlist(list, linecount):
    if linecount:
        list = list[:linecount]
        for x in list:
            print x
    else:
        for x in list:
            print x
    return

def writelist(list, filepath,linecount, filesep,addnewline):
    print '\nwritelist started...'
    wfile = open(filepath, 'w')
    if linecount:
        list = list[:linecount]
    for line in list:
        templine = []
        for elem in line:
            templine.append(str(elem))
        if addnewline:
            templine[len(templine)-1]+='\n'    
        wfile.write(filesep.join(templine))
    wfile.close()
    print '\t wrote '+str(len(list))+' lines'
    print 'writelist ended...'
    return

"""
@param filepath:path to file
@param filesep:delimiter for the columns
@param genpos:gen column position counting from left to right starting at 0
@param gensep:delimiter for multiple genes in one column
@param linecount:set 0 for end of file 
"""
def expand_and_delete(filepatha,filepathb,delimitera, genpos, gensep):
    print '\nexpanding and delete'
    o = open(filepatha,'r')
    n = open (filepathb,'w')
    n.write(o.readline())
    i = 0
    for line in o:
        line = line.split(delimitera)
        line[genpos] = line[genpos].replace('"','')
        if line[genpos].count(gensep) >= 1 :
            temp_dup = line[genpos].split(gensep)
            tempc=[]
            for v in temp_dup:
                tempc.append(v.strip())
                
            temp_nodup = list(sets.Set(tempc))   
            if len(temp_nodup) == 1 :
                line[genpos] = temp_nodup.pop()
                n.write(delimitera.join(line))
            else :
                
                for elem in temp_nodup :
                    line[genpos] = elem
                    #auf grund nebenlaeufigkeit muss diese neuberechnung erfolgen
                    temp = ','.join(line)
                    tempn = temp.split(',')
                    n.write(delimitera.join(tempn))
        else :
            if line[genpos] != '':
                line[genpos]=line[genpos].strip()
                n.write(delimitera.join(line))
        i += 1
        
    print '\twrote lines\n...expand_and_delete done\n'
    
    return

def listtodict(array, genpos):
    print '\nlist to dict started...'
    print '\tlaenge vor dict : ' + str(len(array))
    array.pop(0)
    dict = {}
    i = 0
    for line in array:
        dict.setdefault(line[genpos], []).append(i)
        i += 1
    print '\tlaenge nach dict : ' + str(len(dict))
    print '...list to dict done'
    
    return dict
        
def methy_genexpr_kombinations(metharray, methpos, genarray, genpos, methvaluetoaddlist, genvaluetoaddlist):
    print '\nmethy_genexpr_kombinations started...'
    #resolve genpositions
    methtitlepos = getgenpos(metharray, methvaluetoaddlist)
    genexprpos = getgenpos(genarray, genvaluetoaddlist)
    temparray = []
    title=[]
    for x in methtitlepos:
        title.append(metharray[0][int(x)])
    for y in genexprpos:
        title.append(genarray[0][int(y)])
    metharray = metharray[1:]
    genarray = genarray[1:]
    temparray.append(title)
    
    methdict = listtodict(metharray, methpos)
    gendict = listtodict(genarray, genpos)


    print 'Dictionary gefuellt' 
    nomatch = 0
    while(len(methdict)):        
        #vielleicht hier ander listen benutzen da popitem die liste neu berechnet
        methgen = methdict.popitem()    
        try :
            gengen = gendict.pop(methgen[0])
            for x in methgen[1]:
                for y in gengen:
                    templine = []
                    for m in methtitlepos:
                        templine.append(metharray[int(x)][m]) 
                    for g in genexprpos:
                        templine.append(genarray[int(y)][g])                                                                                   
           
                    #             Gene            TargetID                     ProbeSetID               MChrom                     MMapinfo                    GStart
                    #temparray.append([methgen[0],metharray[int(x)][1],genarray[int(y)][1],metharray[int(x)][2],metharray[int(x)][3],genarray[int(y)][15]+'\n'])
                    temparray.append(templine)     
        except KeyError:
            nomatch += 1
            pass     
    print 'match fail counter: ' + str(nomatch)
    print '...methy_genexpr_kombinations done'

    return temparray
"""
@requires: files must be downward sorted in the specific column
"""
def striparray(array,addnewline):
    i = 0
    while i < len(array):
        array[i] = array[i].replace('\n','')
        i+=1
    if addnewline:
        array[len(array)-1] = array[len(array)-1]+'\n' 
    return array

def appendcolbyfile(filepatha, filepathb, cola, colb, genpos,newfile):
    print '\nappendcol started' 
    filea = open(filepatha,'r')
    fileb = open(filepathb,'r')    
    filenew = open(newfile,'w')
    head = filea.readline().replace('\n','').split('\t')
    #head = head.replace('\t\n','').split('\t')
    headb = fileb.readline().replace('\n','').split('\t')
    head.extend(headb[colb + 1:])
    """
    pos = findvalinarray(head,'\n')
    head = replaceinarray(head, pos, '\n')
    """
    head[len(head)-1] += '\n'
    
    #head = striparray(head,1)
    
    filenew.write('\t'.join(head))  
    d = {}
    arrayb = fileb.readline().replace('\n','').split('\t')
    i=1
    for linea in filea:               
        linea=linea.split('\t')                
        if  d.has_key(linea[cola]): 
            linea='\t'.join(linea).replace('\n','').split('\t')
            linea.extend(d.get(linea[cola]))                    
            linea[len(linea)-1] += '\n'              
            filenew.write('\t'.join(linea))
            i+=1           
        else:
            # dictionary leeren id in arrayb suchen , id in cache[0] schreiben messwerte in cache[1] schreiben
            d.clear()            
            # falls ids in beiden listen sotiert zaehler laufen lassen und indexssuche in arrayb einschraenken
            #while((s+ib)<=leng and linea[cola]>=arrayb[s+ib][colb]):      
            while(linea[cola]>=arrayb[colb] and arrayb[colb] ):            
                if(linea[cola]==arrayb[colb]):   
                    
                    temp = arrayb[colb + 1:]
                    temp ='\t'.join(temp).replace('\n','').split('\t')
                    d.setdefault(linea[cola],temp)  
                    linea='\t'.join(linea).replace('\n','').split('\t')                   
                    linea.extend(d.get(linea[cola])) 
               
                    #linea = striparray(linea,1)
                    linea[len(linea)-1] += '\n'
                      
                    filenew.write('\t'.join(linea))
                    i+=1                                  
                    break 
                arrayb = fileb.readline().split('\t') 
    print str(i)+' lines written to '+str(newfile)    
    filea.close()
    fileb.close()    
    filenew.close()      
    return

def appendcolbyarray(arraya, arrayb, cola, colb, genpos):
    print 'appendcol started'
    temparray = []    
    d = {}
    temparray.append(arraya[0])
    temparray[0].extend(arrayb[0][colb + 1:])
    s=0
    arrayb.pop(0)
    leng = len(arrayb)-1
    for linea in arraya[1:]:
        
        if  d.has_key(linea[cola]):                       
            templine = []
            templine.extend(linea)
            templine.extend(d.get(linea[cola]))
            temparray.append(templine) 
        else:
            # dictionary leeren id in arrayb suchen , id in cache[0] schreiben messwerte in cache[1] schreiben
            d.clear()
            # falls ids in beiden listen sotiert zaehler laufen lassen und indexssuche in arrayb einschraenken
            ib=0
            #print str(linea[cola])+' '+str(arrayb[s+ib][colb])
            while((s+ib)<=leng and linea[cola]>=arrayb[s+ib][colb]):
               
                if(linea[cola]==arrayb[s+ib][colb]):
                                   
                    d.setdefault(linea[cola],arrayb[s+ib][colb + 1:]) 
                    templine = []
                    templine.extend(linea)
                    templine.extend(d.get(linea[cola]))
                    temparray.append(templine)                    
                    ib+=1
                    s += ib
                    break
                ib+=1       
            else:
                s+=ib                
                pass            
           
    return temparray


def sortlist(list,colpos,withhead):
    print '\nsorting list at column '+str(colpos)
    if withhead:
        temp = list.pop(0)
        list = sorted(list, key=lambda list:list[colpos])
        list.insert(0, temp)
    else:
        list = sorted(list, key=lambda list:list[colpos])
    return list

def taketime(f,c,*args):    
    import time 
    start= time.clock()
    for i in range(c):
        f(*args)
    end= time.clock()
    time= (end-start)/c    
    print '\ndurchlaufe: '+str(c) +'\nzeit pro durchlauf: '+str(time) +'\nzeit gesamt: '+str(c*time)
    return time  

def pearsonrandmean(x, y):
    from itertools import imap
    # Assume len(x) == len(y)
    n = len(x)
    sum_x = float(sum(x))
    sum_y = float(sum(y))
    sum_x_sq = sum(map(lambda x: pow(x, 2), x))
    sum_y_sq = sum(map(lambda x: pow(x, 2), y))
    psum = sum(imap(lambda x, y: x * y, x, y))
    num = psum - (sum_x * sum_y/n)
    den = pow((sum_x_sq - pow(sum_x, 2) / n) * (sum_y_sq - pow(sum_y, 2) / n), 0.5)
    if den == 0: return 0
    return num / den

def addcalcpearsonandmean(rpath,delimiter,xpos,ypos,wpath):
    print '\naddcalcpearsonandmean started'
    rfile = open(rpath,'r')
    wfile = open(wpath,'w')
    head=rfile.readline().replace('\n','').split(delimiter)
    head.append('Cor\n')
    wfile.write(delimiter.join(head))    
    for line in rfile:
        line = line.replace('\n','').split(delimiter)
        #print '\n'
        xval=[]
        yval=[]        
        ilist = range(len(xpos))        
        for i in ilist:
            if (line[xpos[i]] != '' and line[xpos[i]] != ' ')and(line[ypos[i]] != '' and line[ypos[i]] != ' '):
                xval.append(float(line[xpos[i]]))
                yval.append(float(line[ypos[i]]))             
            """
            if (line[xpos[i]] != '' and line[xpos[i]] != ' ')and(line[ypos[i]] != '' and line[ypos[i]] != ' '):
                xval.append(float(line[xpos[i]].replace(',','.')))
                yval.append(float(line[ypos[i]].replace(',','.')))  
            """
        #print str(len(xval))+' '+str(xval)
        #print str(len(yval))+' '+str(yval)
        line.append(pearsonrandmean(xval, yval))
        #print corandmean
        line[len(line)-1] = str(line[len(line)-1])+'\n'
        wfile.write(delimiter.join(line))  
    rfile.close()
    wfile.close()
    return 

def replacvalueinearray(array,old,new):
    print '\nreplacing '+str(old)+' with '+str(new)
    head = array.pop(0)
    i=0
    while i < len(array):        
        j = 0
        while j < len(array[i]):
            array[i][j] = array[i][j].replace(old,new)
            j+=1
        i+=1
    array.insert(0,head)
    return array

def sortfilestream(oldfile,writefile,delimiter,sortlist):
    print '\nsorting file as stream'
    file = open(oldfile,'r')
    b = sortlist
    newfile = open(writefile,'w')
  
    try :
        for line in file:            
            line = line.split(delimiter)
            oldline = list(line)
            i=0
            for x in b :
                line[x[0]]=oldline[x[1]]                
                i+=1             
            temp = delimiter.join(line)
            newfile.write(temp)
    except IndexError:
        print 'Row Empty ?! File was written !!!'
        pass
   
    file.close()
    newfile.close()
    print 'done'
    return


def sortgenexpr():
    """
    file = open('CD71GenexpressionsrohdatenNEUGENERIERT.TXT','r')
    head = file.readline().split('\t')
    print len(head)
    file.close()
    sortlist = [[1,getgenpos(head,'0791')],[2,getgenpos(head,'4_0874-08')],[3,getgenpos(head,'0826-08')],[4,getgenpos(head,'5_0890-08')],[5,getgenpos(head,'0834-08')],[6,getgenpos(head,'2_0847-08')],[7,getgenpos(head,'0940-08')],[8,getgenpos(head,'0856-08')],[9,getgenpos(head,'6_0944-08')],[10,getgenpos(head,'3_0861-08')],[11,getgenpos(head,'0952-08')],[12,getgenpos(head,'0980-08')],[13,getgenpos(head,'1062-08')],[14,getgenpos(head,'1079-08')],[15,getgenpos(head,'1307-09')],[16,getgenpos(head,'9_0315-09')],[17,getgenpos(head,'10_0321-09')],[18,getgenpos(head,'11_0327-09')],[19,getgenpos(head,'7_0303-09')],[20,getgenpos(head,'8_0309-09')],[21,getgenpos(head,'0791p-value')],[22,getgenpos(head,'4_0874-08p-value')],[23,getgenpos(head,'0826-08p-value')],[24,getgenpos(head,'5_0890-08p-value')],[25,getgenpos(head,'0834-08p-value')],[26,getgenpos(head,'2_0847-08p-value')],[27,getgenpos(head,'0940-08p-value')],[28,getgenpos(head,'0856-08p-value')],[29,getgenpos(head,'6_0944-08p-value')],[30,getgenpos(head,'3_0861-08p-value')],[31,getgenpos(head,'0952-08p-value')],[32,getgenpos(head,'0980-08p-value')],[33,getgenpos(head,'1062-08p-value')],[34,getgenpos(head,'1079-08p-value')],[35,getgenpos(head,'1307-09p-value')],[36,getgenpos(head,'0903-08p-value')],[37,getgenpos(head,'10_0321-09p-value')],[38,getgenpos(head,'11_0327-09p-value')],[39,getgenpos(head,'7_0303-09p-value')],[40,getgenpos(head,'8_0309-09p-value')],[41,getgenpos(head,'0903-08')],[42,getgenpos(head,'9_0315-09p-value\n')]]        
    print len(sortlist)
    print sortlist
    sortfilestream('CD71GenexpressionsrohdatenNEUGENERIERT.TXT','CD71GenexpressionsrohdatenNEUGENERIERT_sortiert.TXT','\t',sortlist)
    return
    """
    """
    file = open('CD71Methylierungsrohdaten_sortiert.txt','r')
    head = file.readline().split('\t')
    print len(head)
    file.close()
    sortlist = [[14,getgenpos(head,'1062_08.AVG_Beta')],[15,getgenpos(head,'1079_08.AVG_Beta')],[16,getgenpos(head,'1307_09.AVG_Beta')],[17,getgenpos(head,'n0315_09.AVG_Beta')],[18,getgenpos(head,'n0321_09.AVG_Beta')],[19,getgenpos(head,'n0327_09.AVG_Beta')]]        
    print len(sortlist)
    print sortlist
    sortfilestream('CD71Methylierungsrohdaten_sortiert.txt','CD71Methylierungsrohdaten_sortiert_new.txt','\t',sortlist)
    return
    """
    """
    @warning: chromosome col must be like chr01 not chr1 and must be sorted by chromosome and start
    """

                
            
    """
            #colposa = getgenpos('CD71GenexpressionsAnnotationBIG.txt',['Gene Symbol','Chromosome','Start'],'\t')    
            #colposb = getgenpos('refGene.txt',['name2','chrom','txStart','txEnd'],'\t')           
            #print (linea[colposa[1]] == lineb[colposb[1]]) and (int(linea[colposa[2]]) <= int(lineb[colposb[3]]))       
            while linea[colposa[1]] == lineb[colposb[1]] and (int(linea[colposa[2]]) >= int(lineb[colposb[2]])):
                if int(linea[colposa[2]]) <= int(lineb[colposb[3]]):
                    print 'a '+str(a)+' '+linea[colposa[1]]+' '+linea[colposa[2]]+'\nb '+str(b)+' '+lineb[colposb[1]]+' '+lineb[colposb[2]]+' '+lineb[colposb[3]]+'\n'
                    linea[colposa[0]] = lineb[colposb[0]]
                    filew.write(delimitera.join(linea))
                lineb = fileb.readline().split(delimiterb)
                b+=1
            else:
                print 'chrom ungleich'
                print 'a '+str(a)+' '+linea[colposa[1]]+' '+linea[colposa[2]]+'\nb '+str(b)+' '+lineb[colposb[1]]+' '+lineb[colposb[2]]+' '+lineb[colposb[3]]+'\n'
    """       
    filea.close()
    fileb.close()
    filew.close()
    return

def fillemptygenescreatefiles(patha,pathb,chrompos,genepos,repval,delimitera):
    filea = open(patha,'r')
    fileb = open(pathb,'r')
    pos = 0
    heada = filea.readline()
    headb = fileb.readline()
    
    chromlist = []
    for linea in filea:
        linea = linea.split(delimitera)
        if linea[chrompos] not in chromlist:
            chromlist.append(linea[chrompos])
    filenamelist = [] 
    print chromlist 
    filea.seek(0)
    for chr in chromlist:
        filename = 'genexpr'+str(chr)+'.tsv'
        filenamelist.append(filename)
        filew = open(filename,'w')
        filew.write(heada)        
        
        filea.seek(0)
        
        for linea in filea:               
               
            linea=linea.split(delimitera)
            if linea[genepos] == repval and linea[chrompos] == chr:
                filew.write(delimitera.join(linea))
            if linea[chrompos] > chr:
                break            
        filew.close()
        print filename +' closed'
    filea.close()
    fileb.close()

def createfilewithcorgreater(pathr,pathw,cor,sep,pos):
    pos = getgenpos('kombinationenMethrohGenexprohCorMean.tsv','Cor','\t')
    file = open('kombinationenMethrohGenexprohCorMean.tsv','r')
    wfile = open('corgr0.4.tsv','w')
    i = 0
    wfile.write(file.readline())
    c = 0
    for line in file:
        if abs(float(line.split('\t')[pos])) > 0.4:
            wfile.write(line)
            i+=1
        if c%100000 == 0:
            print c
        c+=1        
    print i
    
def addzero(pathr,pathw,delimiter,colpos,valpos=None,val=None):
    filer = open(pathr,'r')
    filew = open(pathw,'w')
    filew.write(filer.readline())
    for line in filer:
        line = line.split(delimiter)
        col = line[colpos]
        
        if len(col) == 4:
            if(line[colpos] != 'chrX' and line[colpos] != 'chrY'):
                line[colpos] = line[colpos][:valpos]+str(val)+line[colpos][valpos:]
        if len(col) in [1,2]:
            if len(col) == 1: 
                if col in ['X','Y']:
                    line[colpos] = 'chr'+col
                else:                    
                    line[colpos] = 'chr0'+col
            else:
                line[colpos] = 'chr'+col
    
        filew.write(delimiter.join(line))
    return                


def calcmean(line,pos):
    sum = 0    
    failcount=0
    for i in pos:
        try:
            #sum+=float(line[i].replace(',','.'))
            sum+=float(line[i])
        except ValueError:
            failcount+=1
    
    try:
        return sum/(len(pos)-failcount)
    except ZeroDivisionError:
        print 'division by zero mean set to 0'
        return 0
    

def addcalcmean(rpath,wpath,delimiter,m1pos,m2pos,m3pos,m4pos):
    print '\naddcalcmean started...'
    rfile = open(rpath,'r')
    wfile = open(wpath,'w')
    head=rfile.readline().replace('\n','').split(delimiter)
    nold = len(head) - 1
    head.extend(['MeanBetaMDS','MeanBetaHealthy','DiffBeta','MeanExpMDS','MeanExpHealthy','DiffExp\n'])
    nnew = len(head) - 1
    wfile.write('\t'.join(head))
    
    for line in rfile:
        line = line.split(delimiter)
        line[nold] = line[nold].replace('\n','')
        m1=calcmean(line,m1pos)
        m2=calcmean(line,m2pos)
        m3=calcmean(line,m3pos)
        m4=calcmean(line,m4pos)
        line.extend([str(m1),str(m2),str(m1-m2),str(m3),str(m4),str(m3-m4)])              
        line = delimiter.join(line)  
        wfile.write(line+'\n')
    
    rfile.close()
    wfile.close()
    print 'addcalcmean ended...'
    return

def addcalcmean2(rpath,wpath,delimiter,mpos,coltitle):
    print '\naddcalcmean started...'
    rfile = open(rpath,'r')
    wfile = open(wpath,'w')
    head=rfile.readline().replace('\n','').split(delimiter)
    nold = len(head) - 1
    head.extend(coltitle)
    nnew = len(head) - 1
    #new
    head[nnew] +='\n'
    wfile.write('\t'.join(head))
    
    for line in rfile:
        line = line.split(delimiter)
        line[nold] = line[nold].replace('\n','')
        for m in mpos:
            mcalc = calcmean(line,m)
            line.append(str(mcalc))                 
        line = delimiter.join(line)  
        wfile.write(line+'\n')
    
    rfile.close()
    wfile.close()
    print 'addcalcmean ended...'
    return

def addcalcdiff(rpath,wpath,delimiter,mpos,coltitle):
    print '\naddcalcdiff started...'
    rfile = open(rpath,'r')
    wfile = open(wpath,'w')
    head=rfile.readline().replace('\n','').split(delimiter)
    nold = len(head) - 1
    head.extend(coltitle)
    nnew = len(head) - 1
    wfile.write('\t'.join(head)+'\n')
    
    for line in rfile:
        line = line.split(delimiter)
        line[nold] = line[nold].replace('\n','')
        for m in mpos:
            line.append(str(float(line[m[0]])-float(line[m[1]])))           
        line = delimiter.join(line)  
        wfile.write(line+'\n')
    
    rfile.close()
    wfile.close()
    print 'addcalcdiff ended...'
    return

def countvalue(path,delimiter,value,pos):
    f = open(path,'r')
    i = 0
    for line in f:
        line = line.split(delimiter)
        if line[pos] == value:
            i=i+1
    f.close()
    return i

def fillemptygenes(patha,pathb,pathw,colposa,colposb,repval,delimitera,delimiterb):
    filea = open(patha,'r')
    fileb = open(pathb,'r')
    filew = open(pathw,'w')
    fileb.readline()
    filew.write(filea.readline())
    
    lineb=fileb.readline().split(delimiterb)
    a=2
    wrote = []
    pos = [0]
    newstart = 0
    for linea in filea:        
        wrote = []      
        linea = linea.split(delimitera)
        print 'a '+str(a)+' '+linea[colposa[1]]+' '+linea[colposa[2]]+'\nb '+str(int(pos[0]))+' '+lineb[colposb[1]]+' '+lineb[colposb[2]]+' '+lineb[colposb[3]]+'\n'
        #colposa = getgenpos('2014/sampleRohdatenGENE-CORE(65)linear.RMA-GENE-CORE-Group2log2SORT.tsv',['Gene Symbol','Chromosome','Start'],'\t')
        #colposb = getgenpos('2014/refGeneSorted.tsv',['name2','chrom','txStart','txEnd'],'\t')

        #if (a > 80):break
        if linea[colposa[0]] == repval:    
            # chr a gleich chr b
            if newstart:
                fileb.seek(newstart)
                #print 'took newstart '+str(newstart)
                lineb = fileb.readline().split(delimiterb)
                newstart = 0           
            if linea[colposa[1]] == lineb[colposb[1]]:
                #print 'a '+str(a)+' '+linea[colposa[1]]+' '+linea[colposa[2]]+'\nb '+str(int(pos[0]))+' '+lineb[colposb[1]]+' '+lineb[colposb[2]]+' '+lineb[colposb[3]]+'\n'

                # so lange start > txstart
                while int(linea[colposa[2]]) >= int(lineb[colposb[2]]):
                    # falls chr a > chr b break , nextline in b               
                    if linea[colposa[1]] < lineb[colposb[1]]:
                        break
                    #gen schon geschrieben ansonsten redundante eintraege
                    if lineb[colposb[0]] in wrote:
                        pass
                        #print 'inwrote !!!'
                        #print '\ta '+str(a)+' '+linea[colposa[1]]+' '+linea[colposa[2]]+'\n\tb '+str(int(pos[0]))+' '+lineb[colposb[1]]+' '+lineb[colposb[2]]+' '+lineb[colposb[3]]
                        
                    if int(linea[colposa[2]]) <= int(lineb[colposb[3]]) and (lineb[colposb[0]] not in wrote):
                        newstart = int(pos[0])
                        wrote.append(lineb[colposb[0]])
                        #print 'a '+str(a)+' '+linea[colposa[1]]+' '+linea[colposa[2]]+'\nb '+str(int(pos[0]))+' '+lineb[colposb[1]]+' '+lineb[colposb[2]]+' '+lineb[colposb[3]]
                        #print wrote
                        linea[colposa[0]] = lineb[colposb[0]]
                        
                        #print ''
                    lineb = fileb.readline().split(delimiterb)
                    if not wrote:
                        pos.append(int(fileb.tell()))
                        if len(pos) > 2:pos.pop(0)
                        #print pos
       
            # chr a < chr b weiter in a glaube das wird nie erreicht
            elif linea[colposa[1]] < lineb[colposb[1]]:
                print 'a '+str(a)+' '+linea[colposa[1]]+' '+linea[colposa[2]]+'\nb '+str(int(pos[0]))+' '+lineb[colposb[1]]+' '+lineb[colposb[2]]+' '+lineb[colposb[3]]
                        
                a+=1
                continue
            # chr a > chr b weiter in b glaube das wird nie erreicht
            elif linea[colposa[1]] > lineb[colposb[1]]:
                lineb = fileb.readline().split(delimiterb)
                pos.append(int(fileb.tell()))
                if len(pos) > 2:pos.pop(0)
            
        filew.write(delimitera.join(linea))
        a+=1


def runit():
    #os.system('SortFile.py refGene.tsv refGene.tsv -k line.split("\t")[2]')
    #addzero('CD71GenexpressionsAnnotationBIG2.txt', 'CD71GenexpressionsAnnotationBIG2zero.txt', '\t', 12, 3, 0)
    #colposa = getgenpos('CD71GenexpressionsAnnotationBIG2zero.txt',['Gene Symbol','Chromosome','Start'],'\t')    
    #colposb = getgenpos('refGeneSorted.tsv',['name2','chrom','txStart','txEnd'],'\t')
    #fillemptygenes('CD71GenexpressionsAnnotationBIG2zero.txt','refGeneSorted.tsv', 'CD71GenexpressionsAnnotationBIG2zero_addedfinal.txt',colposa,colposb,'---', '\t','\t')          
    
      

    #addzero('refGene.txt','refGene.tsv','\t',getgenpos('refGene.txt', 'chrom','\t'),3,'0')
    #addzero('HumanMethylation450.csv','HumanMethylation450_chr.csv',',', getgenpos('HumanMethylation450.csv','CHR',','))  
    #a = readfile('HumanMethylation450_chr.csv', ',', 0)
    #writelist(a,'HumanMethylation450_chr.tsv', 0, '\t',0)
    
    meth = readfile('HumanMethylation450.csv',',',0)
    #genexpr = readfile('CD71Genexpressionsannotation_noemptygenesID.tsv','\t',0)  
    genexpr = readfile('CD71GenexpressionsAnnotationBIG2zero_addedfinal.txt','\t',0)
     
    meth = stripper(meth,getgenpos(meth,'UCSC_RefGene_Name'))  
    meth = expand_and_delete(meth,getgenpos(meth,'UCSC_RefGene_Name'),';')  
    
    #genexpr = readfile('CD71GenexpressionsAnnotationBIG.txt','\t',0)   
    #refGene = readfile('refGene.txt', '\t',0)
    # total falsch meth = fillemptygenes(genexpr,refGene,getgenpos(genexpr,['Gene Symbol','Chromosome','Start']),getgenpos(refGene,['name2','chrom','txStart','txEnd']))  
    
    kombi = methy_genexpr_kombinations(meth,getgenpos(meth,'UCSC_RefGene_Name'),genexpr,getgenpos(genexpr,'Gene Symbol'),['UCSC_RefGene_Name','CHR','IlmnID','MAPINFO'],['Probe Set ID','Start','Level'])    
    del meth
    del genexpr
    writelist(kombi,'kombinationen.tsv',0, '\t', 1)

    #old way to sort
    #kombi = sortlist(kombi,getgenpos(kombi, 'IlmnID'),1)
    #writelist(kombi,'kombinationen.tsv',0, '\t', 1)
    del kombi
    SortFile.main('kombinationen.tsv','kombinationen.tsv','line.split("\t")[2]',64000)
    
    
    #readfile methraw
    #methraw = replacvalueinearray(methraw,',','.')
    #methraw = sortlist(methraw,getgenpos(methraw, 'TargetID'),1)
    #write methraw
    #methraw = readfile('CD71Methyrohdaten_TargetIDsorted.txt', '\t', 0)
    #kombi = appendcolbyarray(kombi,methraw, getgenpos(kombi, 'IlmnID'), getgenpos(methraw, 'TargetID'), getgenpos(kombi, 'UCSC_RefGene_Name'))    
    appendcolbyfile('kombinationen.tsv','CD71Methyrohdaten_TargetIDsorted.txt', getgenpos('kombinationen.tsv', 'IlmnID','\t'), getgenpos('CD71Methyrohdaten_TargetIDsorted.txt', 'TargetID','\t'), getgenpos('kombinationen.tsv', 'UCSC_RefGene_Name','\t'),'kombinationenMethroh.tsv')
    # hier wird die sortierte benutzt evlt vorher sotieren
   
    #genexprraw = readfile('CD71GenexpressionsrohdatenNEUGENERIERT_sortiert_new.TXT', '\t', 0)
    #genexprraw = sortlist(genexprraw,getgenpos(genexprraw, 'Probe Set ID'),1)  
    
    pos = getgenpos("kombinationenMethroh.tsv","Probe Set ID","\t")
    print pos
    SortFile.main('kombinationenMethroh.tsv','kombinationenMethroh.tsv','line.split("\t")[4]',256000)  
    #kombi = readfile('kombinationenMethroh.tsv','\t', 0)
    #kombi = sortlist(kombi,getgenpos(kombi, 'Probe Set ID'),1)  
    #writelist(kombi,'kombinationenMethroh.tsv',0, '\t',0)

    appendcolbyfile('kombinationenMethroh.tsv', 'CD71Genexpr_ProbeSetIDsorted.tsv', getgenpos('kombinationenMethroh.tsv', 'Probe Set ID','\t'), getgenpos('CD71Genexpr_ProbeSetIDsorted.tsv', 'Probe Set ID','\t'), getgenpos('kombinationenMethroh.tsv', 'UCSC_RefGene_Name','\t'),'kombinationenMethrohGenexproh.tsv')        
    
    #kombi = deleteinstring(kombi,getgenpos(kombi,'n0309_09.AVG_Beta\n'),'\n')    
    #kombi = deleteinstring(kombi,getgenpos(kombi,'9_0315-09p-value\n'),'\n')
    
    xpos = getgenpos('kombinationenMethrohGenexproh.tsv',['0791_08.AVG_Beta','0874_08.AVG_Beta','0826_08.AVG_Beta','0890_08.AVG_Beta','0834_08.AVG_Beta','0847_08.AVG_Beta','0940_08.AVG_Beta','0865_08.AVG_Beta','0944_08.AVG_Beta','0861_08.AVG_Beta','0952_08.AVG_Beta','0980_08.AVG_Beta','1062_08.AVG_Beta','1079_08.AVG_Beta','1307_09.AVG_Beta','n0315_09.AVG_Beta','n0321_09.AVG_Beta','n0327_09.AVG_Beta', 'n0303_09.AVG_Beta', 'n0309_09.AVG_Beta'],'\t')
    ypos = getgenpos('kombinationenMethrohGenexproh.tsv',['0791','4_0874-08','0826-08','5_0890-08','0834-08','2_0847-08','0940-08','0856-08','6_0944-08','3_0861-08','0952-08','0980-08','1062-08','1079-08','1307-09','9_0315-09','10_0321-09','11_0327-09','7_0303-09','8_0309-09'],'\t')
    
    addcalcpearsonandmean('kombinationenMethrohGenexproh.tsv','\t', xpos, ypos,'kombinationenMethrohGenexprohCor.tsv')
    
    m1pos=getgenpos('kombinationenMethrohGenexprohCor.tsv',['0940_08.AVG_Beta','0865_08.AVG_Beta','0944_08.AVG_Beta','0861_08.AVG_Beta','0952_08.AVG_Beta','0980_08.AVG_Beta','1062_08.AVG_Beta','1079_08.AVG_Beta','1307_09.AVG_Beta'],'\t')
    m2pos=getgenpos('kombinationenMethrohGenexprohCor.tsv',['n0315_09.AVG_Beta','n0321_09.AVG_Beta','n0327_09.AVG_Beta','n0303_09.AVG_Beta','n0309_09.AVG_Beta'],'\t')
    m3pos=getgenpos('kombinationenMethrohGenexprohCor.tsv',['0940-08','0856-08','6_0944-08','3_0861-08','0952-08','0980-08','1062-08','1079-08','1307-09'],'\t')
    m4pos=getgenpos('kombinationenMethrohGenexprohCor.tsv',['9_0315-09','10_0321-09','11_0327-09','7_0303-09','8_0309-09'],'\t')
        
    addcalcmean('kombinationenMethrohGenexprohCor.tsv','kombinationenMethrohGenexprohCorMean.tsv','\t', m1pos, m2pos, m3pos, m4pos)
    

    #7528891 kombinationen bei neuen files

    return

def oldcalc():
    """
    m1pos=getgenpos('kombinationenMethrohGenexprohCor.tsv',['0940_08.AVG_Beta','0865_08.AVG_Beta','0944_08.AVG_Beta','0861_08.AVG_Beta','0952_08.AVG_Beta','0980_08.AVG_Beta','1062_08.AVG_Beta','1079_08.AVG_Beta','1307_09.AVG_Beta'],'\t')
    m2pos=getgenpos('kombinationenMethrohGenexprohCor.tsv',['n0315_09.AVG_Beta','n0321_09.AVG_Beta','n0327_09.AVG_Beta','n0303_09.AVG_Beta','n0309_09.AVG_Beta'],'\t')
    m3pos=getgenpos('kombinationenMethrohGenexprohCor.tsv',['0940-08','0856-08','6_0944-08','3_0861-08','0952-08','0980-08','1062-08','1079-08','1307-09'],'\t')
    m4pos=getgenpos('kombinationenMethrohGenexprohCor.tsv',['9_0315-09','10_0321-09','11_0327-09','7_0303-09','8_0309-09'],'\t')
    headtitle = ['MeanBetaMDS','MeanBetaHealthy','MeanExpMDS','MeanExpHealthy']
    #headtitle = ['MeanBetaMDS','MeanBetaHealthy','MeanExpMDS','MeanExpHealthy\n']
    addcalcmean2('kombinationenMethrohGenexprohCor.tsv','kombinationenMethrohGenexprohCorMean.tsv','\t',[m1pos, m2pos, m3pos, m4pos],headtitle)

    m1pos=getgenpos('kombinationenMethrohGenexprohCorMean.tsv','MeanBetaMDS','\t')
    m2pos=getgenpos('kombinationenMethrohGenexprohCorMean.tsv','MeanBetaHealthy','\t')
    m3pos=getgenpos('kombinationenMethrohGenexprohCorMean.tsv','MeanExpMDS','\t')
    m4pos=getgenpos('kombinationenMethrohGenexprohCorMean.tsv','MeanExpHealthy\n','\t')
    headtitle = ['DiffMeanBetaMDS_MeanBetaHealthy','DiffMeanExpMDS_MeanExpHealthy']
    addcalcdiff('kombinationenMethrohGenexprohCorMean.tsv','kombinationenMethrohGenexprohCorMeanDiff.tsv','\t',[[m1pos, m2pos],[m3pos, m4pos]],headtitle)
    """
    #appendcolbyfile('test/GenExprAnno.tsv', 'test/GenExprRoh.tsv', 0, 0, 1, 'test/Ergebnis.tsv')
    #appendcolbyfile('test/CD71GenexpressionsAnnotationBIG2zero_addedfinal.txt', 'test/CD71Genexpr_ProbeSetIDsorted.tsv', 0, 0, 1, 'test/CD71GenexprAnnoNRoh.tsv')
    #printfile('test/CD71GenexpressionsAnnotationBIG2zero_addedfinal.txt', 10, '\t')
    #SortFile.main('test/CD71GenexprAnnoNRoh_small.tsv','test/CD71GenexprAnnoNRoh_small.tsv','str(line.split("\t")[1])',128000)
    #print getgenpos('corgr0.4.tsv','Cor','\t')
    #SortFile.main('corgr0.4.tsv','corgr0.4sort.tsv','abs(float(line.split("\t")[68]))',128000)

def newcalc():
    #merge new meth microarray files with old microarray files - do it bot ways to ensure max overlap
    #appendcolbyfile('2014/CD71Methyrohdaten_TargetIDsortedALT.txt','2014/131203_Meth125_Mossner_SampleMethylationProfile_all_AnnoSort.tsv', 1, 0, 3,'2014/methmatch1.tsv')
    #appendcolbyfile('2014/131203_Meth125_Mossner_SampleMethylationProfile_all_AnnoSort.tsv','2014/CD71Methyrohdaten_TargetIDsortedALT.txt', 0, 1, 3,'2014/methmatch2.tsv')
    return

def calcLogLinear(filea,fileb,startcol):
    print 'calcLogLinear started'
    f = open(filea,'r')
    fw = open(fileb,'w')
    fw.write(f.readline())
    for line in f:
        line = line.strip().split('\t')
        value = line[startcol:]
        
        i=0
        for i in range(len(value)):
            value[i]=float(value[i])
        
        value = np.power(2,value)
        
        erg = []
        erg.extend(line[:startcol])
        
        for v in value:
            erg.append(str(round(v,4)))

        erg = '\t'.join(erg)+'\n'
        
        fw.write(erg)
       
    f.close()
    fw.close()
    print 'calcLogLinear ended'
    
    

def calccornew():
    zu = open('2014/zuordnungFuerKorrelation.txt','r')
    zu.readline()
    zu.readline()
    exp = zu.readline().strip().split('\t')
    meth = zu.readline().strip().split('\t')
    zu.close()
    o = open('2014/Cor_CD71MethExpCore.tsv','r')
    n = open('2014/Cor_CD71MethExpCore_R.tsv','w')
    n.write(o.readline().strip()+'\tCor\tCor_PValue\n')
    for line in o:
        try:
            line = line.strip().split('\t')
            x = []
            y=[]
            for v in exp:
                x.append(float(line[int(v)]))
            for v in meth:
                y.append(float(line[int(v)]))

            temp = stats.pearsonr(x, y)
            line.append(str(temp[0]))
            line.append(str(temp[1]))
            n.write('\t'.join(line)+'\n')
        except ValueError:
            continue

    o.close()
    n.close()

def fillemptynew(patha,pathb,pathw,pathwe,gencola,gencolb,starta,startb,repval):
    print 'fillemptynew started'
    """
    replace unannotated gene entries like --- in path a. search in file b for matching transcript position
    """
    fa = open(patha,'r')
    fr = open(pathb,'r')    
    fw = open(pathw,'w')
    fwe = open(pathwe,'w')
    
    
    anotitle = fa.readline()
    fw.write(anotitle)
    
    d=[]
    reftitle = fr.readline().split('\t')
    refgenecol = reftitle.index('name2')
    refstart = reftitle.index('txStart')
    refstop = reftitle.index('txEnd')
    refchrom = reftitle.index('chrom')
    
    chromdict = {}
    

    for line in fr:
       
        line = line.split('\t')
        temp = []
        temp.append(line[refstart])
        temp.append(line[refstop])
        temp.append(line[refgenecol])  
        temp.append(line[refchrom])      
        
        #-1 because array chr1 starts at array pos 0
        if line[refchrom] == 'chrX' or line[refchrom] == 'chrY':
            if line[refchrom] == 'chrX':
                chromdict.setdefault(line[refchrom],23-1)
            else:
                chromdict.setdefault(line[refchrom],24-1)
            
        else:
            chromdict.setdefault(line[refchrom],int(line[refchrom][3:])-1)
            
        d.append(temp)
            
    
    k = sorted(d,key=lambda x:int(x[1]))
    k = sorted(d,key=lambda x:int(x[0]))
    #k0: start k1: stop k2: gen k3: chrom
    print len(k)
    print len(chromdict)
    print chromdict
    
    chromseplist = []
    
    itter = []
    #itter from chr1 to chr 22 + chrX + chrY
    for i in range(1,len(chromdict)-1):
        itter.append('chr'+str(i))
    itter.append('chrX')
    itter.append('chrY')
    
    print itter
    
    # filter k by chromosom create new list for each chromosom. append new array with chromosoms to chomlist
    for s in itter:    
        temp = filter(lambda v:  v[3] == s, k) 
        chromseplist.append(temp)
        #just a check
        #fw.write(str(temp)+'\n') 
    
    anotitle = anotitle.split('\t')
    anogenecol = anotitle.index('Gene Symbol')
    anostart = anotitle.index('Start')
    anostop = anotitle.index('Stop')
    anochrom = anotitle.index('Chromosome')
    
    
    i=0
    for anoline in fa:
        try:
            anoline = anoline.split('\t')
            #search the appropriate chromosome list
            found = False
            for refline in chromseplist[chromdict.setdefault(anoline[anochrom])]:            
                if int(refline[0]) <= int(anoline[anostart]) and int(refline[1]) >= int(anoline[anostop]) and refline[3].strip() == anoline[anochrom].strip():
                    anoline[anogenecol] = refline[2]
                    fw.write('\t'.join(anoline)) 
                    found = True
                    
            if not found:
                fwe.write('\t'.join(anoline))                  
                    
            i=i+1
            if i%10000 == 0:
                print str(i/10000)+str('% processed')
        except TypeError:
            fwe.write('EXCEPTION\n')
            fwe.write('\t'.join(anoline))
   
     
    fa.close()
    fw.close()
    fr.close()
    """
    #call
    #FillemptygenesNew
    filea = '2014All/Edit/AdditionalRefGeneAnnotation/sampleRohdatenEXON-ALL(65)ohneYundControlsQlucoreRMA-EXON-ALL-DABGGroup9LinearNoANNO.tsv'
    fileb = 'refgene/refGene.tsv'
    filec = '2014All/Edit/AdditionalRefGeneAnnotation/sampleRohdatenEXON-ALL(65)ohneYundControlsQlucoreRMA-EXON-ALL-DABGGroup9LinearAddANNO.tsv'
    filed = '2014All/Edit/AdditionalRefGeneAnnotation/sampleRohdatenEXON-ALL(65)ohneYundControlsQlucoreRMA-EXON-ALL-DABGGroup9LinearAddANNOErrors.tsv'
    #filec = '2014All/Edit/AdditionalRefGeneAnnotation/sampleRohdatenEXON-ALL(65)ohneYundControlsQlucoreRMA-EXON-ALL-DABGGroup9LinearGenesSort.tsv'  
    fillemptynew(filea,fileb,filec,filed,'','', '','','')
    
    #call to eliminate duplicated gen entries
    
    filea = '2014All/Edit/AdditionalRefGeneAnnotation/sampleRohdatenEXON-ALL(65)ohneYundControlsQlucoreRMA-EXON-ALL-DABGGroup9LinearAddANNO.tsv'    
    
    filew = '2014All/Edit/AdditionalRefGeneAnnotation/sampleRohdatenEXON-ALL(65)ohneYundControlsQlucoreRMA-EXON-ALL-DABGGroup9LinearAddANNONoMulti.tsv'
    #filec = '2014All/Edit/AdditionalRefGeneAnnotation/sampleRohdatenEXON-ALL(65)ohneYundControlsQlucoreRMA-EXON-ALL-DABGGroup9LinearGenesSort.tsv'  
    delmultiplelines(filea,filew)
    
    """
    
def delmultiplelines(filepatha,filepathb):
    f = open(filepatha,'r')
    fw = open(filepathb,'w')
    fw.write(f.readline())
    
    #offset = f.tell()
    dict = {}
    last = f.readline().split('\t')
    
    i=0
    while True:
        try:
            line = f.readline()        
            line = line.split('\t')     
                   
            if line[1] == last[1]:
                temp = []
                temp.extend(line)
                dict.setdefault(line[0],temp)
                
            else:
                if len(dict) == 1:                      
                    fw.write('\t'.join(last))
                    dict = {}
                    temp = []
                    temp.extend(line)
                    dict.setdefault(line[0],temp)
                
                else:
                    temp = []
                    for v in dict:
                        temp.append(dict.setdefault(v))
                    #  k = sorted(d,key=lambda x:int(x[1]))
                    #templist = sorted(temp, key=lambda x:int(x[0]))
                    for v in temp:
                        fw.write('\t'.join(v))
                        
                    dict = {}
                    temp = []
                    temp.extend(line)
                    dict.setdefault(line[0],temp)
            
            #offset.tell()            
            last = line
        except IndexError:
            temp = []
            for v in dict:
                temp.append(dict.setdefault(v))
            for v in temp:
                fw.write('\t'.join(v))
            break
            
    f.close()
    fw.close()
    print 'done'

def calcSplicingIndex(filepatha,filepathb,filepathc):
    """
    @param filepatha: gene core file with genewide expression for each patient. NO PROBESETS INCLUDED
    @param filepathb: probesetexpression file with patients
    @param filepathc: write new file to filepathc
    """
    fc = open(filepatha,'r')
    fp = open(filepathb,'r')
    fw = open(filepathc,'w')
    fwe = open(filepathc+'errors','w')
    genecoretitles = fc.readline().strip().split('\t')
    probesettitles = fp.readline().strip().split('\t')
    
    genexprdict = {}
    for line in fc:
        line = line.strip().split('\t')
        genexprdict.setdefault(line[0],line)
    
    
    mapper2 = File.Mapper('maps\SplicingTranscriptExpressionMapping.csv',':');
    
    valtitles = mapper2.extractTitles()
    
    posmapper = {}
    
    for p in probesettitles:
        if p in valtitles:
            #key:new pos val:old pos
            posmapper.setdefault(probesettitles.index(p),genecoretitles.index(p))
     
    poslist = []     
    for p in probesettitles:
        if p in valtitles:
            temp = []
            temp.append(probesettitles.index(p))
            temp.append(genecoretitles.index(p))
            poslist.append(temp)
            
    print 'posmapper %s'%posmapper
    print len(posmapper)
    
    print 'poslist %s'%poslist
    print len(poslist)
    
    print 'valtitles %s'%valtitles
    
    print 'groupmaps'
    groupmaps = mapper2.extractGroupMapping()
    print groupmaps
    
    newtitles =[]
    
    for p in probesettitles:     
        if p in valtitles:
            newtitles.append(p+'_NI')
            
    print 'newtitles %s'%newtitles
    
    newprobesettitles = []
    newprobesettitles.extend(probesettitles)
    newprobesettitles.extend(newtitles)
    
    print 'newprobesettitles %s'%newprobesettitles
    
    dict = mapper2.extractGroupMapping()
    pos = mapper2.findPositions(valtitles,mapper2.mapperdict)
    
    print 'dict\t%s'%dict
    
        
    newdict = {}
    #dict = mapper2.extractGroupMapping()
    #pos = mapper2.findPositions(valtitles,mapper2.mapperdict)
    
    for v in dict:
        group = dict.setdefault(v)
        for g in group:
            #hexenwerk
            newdict.setdefault(v,[]).append(g+'_NI')
            
    print 'newdict %s' %newdict
    
    newpos = mapper2.findPositions(newprobesettitles,newdict)
    #print 'pos\t%s'%pos    
    for i in range(len(newpos[0])):newpos[0][i]+='_NI'
    #print 'newpos\t%s' %newpos
    for v in pos:print v
    for v in newpos:print v
    
    newprobesettitles.extend(newpos[0])    
    SIgrouptitles = [
'SI_healthyNI_VS_lowNI',
'SI_healthyNI_VS_highNI',
'SI_lowNI_VS_highNI',
'SI_preLenNI_VS_postLenNI',
'SI_pre5AzaNI_VS_post5AzaNI',
'SI_pre5AzaNI_VS_preLenNI',
'SI_healthyNI_VS_preLenNI',
'SI_healthyNI_VS_pre5AzaNI']
    
    newprobesettitles.extend(SIgrouptitles)
    
    SIgroups = []
    #group positions for all low patients in one array
    templow = []
    templow.append(newprobesettitles.index('preLen_NI'))
    templow.append(newprobesettitles.index('postLen_NI'))
    SIgroups.append([[newprobesettitles.index('healthy_NI')],templow])
    
    temphigh = []
    temphigh.append(newprobesettitles.index('pre5Aza_NI'))
    temphigh.append(newprobesettitles.index('post5Aza_NI'))    
    SIgroups.append([[newprobesettitles.index('healthy_NI')],temphigh])
    
    SIgroups.append([templow,temphigh])
    
    SIgroups.append([[newprobesettitles.index('preLen_NI')],[newprobesettitles.index('postLen_NI')]])
    
    SIgroups.append([[newprobesettitles.index('pre5Aza_NI')],[newprobesettitles.index('post5Aza_NI')]])
    
    SIgroups.append([[newprobesettitles.index('pre5Aza_NI')],[newprobesettitles.index('preLen_NI')]])
    
    SIgroups.append([[newprobesettitles.index('healthy_NI')],[newprobesettitles.index('preLen_NI')]])
    
    SIgroups.append([[newprobesettitles.index('healthy_NI')],[newprobesettitles.index('pre5Aza_NI')]])


    fw.write('\t'.join(newprobesettitles)+'\n')
    fwe.write('\t'.join(newprobesettitles)+'\n')
    
    i=0
    genecorline = None
    for line in fp:
        line = line.strip().split('\t')        
        
        if genecorline == None or genecorline[0] != line[1]:    
            genecorline = genexprdict.setdefault(line[1])    
            if not genecorline:
                #print 'no match for '+str(line[1])+' skipping to next Gene'
                fwe.write('\t'.join(line)+'\n')
                continue                
        """
        calculate NI for every Patient
        """
        for p in poslist:
                #probeset / gene wide expression line[p[0]]: actual probeset intensity for one patient genecorline[p[1]]: gene wide expression for same patient from genecorefile
                line.append(str(round(float(line[p[0]])/float(genecorline[p[1]]),4)))
        
        """
        calculate group NIs
        """
        for g in newpos[1:]:
            temp=[]
            for v in g:
                temp.append(float(line[v]))
            line.append(str(round(np.mean(temp),4)))
        
        """
        Calculate group SI
    
        healthy_NI_VS_MDSlow+int-1_NI
        healthy_NI_VS_MDSint-2+high_NI
        MDSlow+int-1_NI_VS_MDSint-2+high_NI
        MDSint-2+high_NI_VS_post-5-Aza_NI
        MDSlow+int-1_NI_VS_post-Len_NI
        """
        
        #for low and high there are multiple column indexes
        for g in SIgroups:
            tempg0 = 0
            tempg1 = 0
            for v in g[0]:
                tempg0+=float(line[v])
            tempg0 = tempg0 / len(g[0])
            
            for v in g[1]:
                tempg1+=float(line[v])
            tempg1 = tempg1 / len(g[1])
            
            line.append(str(round(tempg0/tempg1,4)))
        
        
        fw.write('\t'.join(line)+'\n')    
            
        i=i+1
        
        if i%10000 == 0:
            print str(i/10000)+str('% processed')

    fc.close()
    fp.close()
    fw.close()
    fwe.close()
    
if __name__ == '__main__':
    #taketime(runit, 1)
    #filea = '2014/EXON-EXTENDED(65)linear.RMA-EXON-EXTENDED-DABG-Group5log2AddAnnoSort_FINAL.tsv'
    #fileb = '2014/EXON-EXTENDED(65)linear.RMA-EXON-EXTENDED-DABG-Group5log2AddAnnoSort_FINAL2.tsv'
    #SortFile.main(filea,fileb,'str(line.split("\t")[1]),int(line.split("\t")[7]),',128000)
    
    filec = '2014All/Edit/sampleRohdatenEXON-ALL(65)ohneYundControlsQlucoreRMA-EXON-ALL-DABGGroup9LinearAddAnnoFINALSorted_new.tsv'
    
    fileb = '2014All/Edit/sampleRohdatenEXON-ALL(65)ohneYundControlsQlucoreRMA-EXON-ALL-DABGGroup9LinearAddAnnoFINALSorted.tsv'
    #calcLogLinear(filea,fileb,9)
    
    SortFile.main(fileb,filec,'str(line.split("\t")[1]),int(line.split("\t")[7])',128000)
    """
    filea = '2014All/Edit/GENE-CORE(65)linear.RMA-GENE-CORE-Group2SortLinear.tsv'
    fileb = '2014All/Edit/sampleRohdatenEXON-ALL(65)ohneYundControlsQlucoreRMA-EXON-ALL-DABGGroup9LinearAddAnnoFINALSorted.tsv'    
    filec = '2014All/Edit/sampleRohdatenEXON-ALL(65)ohneYundControlsQlucoreRMA-EXON-ALL-DABGGroup9LinearAddAnnoFINALSortedSI.tsv' 
    calcSplicingIndex(filea,fileb,filec)
    """
    #SortFile.main(filea,filec,'str(line.split("\t")[1]),int(line.split("\t")[0])',128000)

    

   