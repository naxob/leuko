import sets
'''
Created on 28.08.2012
@author: Axel Wilbertz
'''
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

def printfile(filepath, linecount):
    file = open(filepath, 'r')
    for line in range(linecount):
        print file.readline()
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
def expand_and_delete(flist, genpos, gensep):
    print '\nexpanding and delete'
    temparray = []
    i = 0
    for line in flist:
        
        if line[genpos].count(gensep) >= 1 :
            temp_dup = line[genpos].split(gensep)
            temp_nodup = list(sets.Set(temp_dup))   
            if len(temp_nodup) == 1 :
                line[genpos] = temp_nodup.pop()
                temparray.append(line)
            else :
                
                for elem in temp_nodup :
                    line[genpos] = elem
                    #auf grund nebenlaeufigkeit muss diese neuberechnung erfolgen
                    temp = ','.join(line)
                    tempn = temp.split(',')
                    temparray.append(tempn)
        else :
            temparray.append(line)   
        i += 1
        
    print '\twrote ' + str(len(temparray)) + ' lines\n...expand_and_delete done\n'
    
    return temparray

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

def appendcolbyfile(filepatha, filepathb, cola, colb, genpos,newfile,*debug):
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
            while(linea[cola]>=arrayb[colb]):               
                if(linea[cola]==arrayb[colb]):   
                    temp = arrayb[colb + 1:]
                    temp ='\t'.join(temp).replace('\n','').split('\t')
                    d.setdefault(linea[cola],temp)  
                    linea='\t'.join(linea).replace('\n','').split('\t')                   
                    linea.extend(d.get(linea[cola])) 
               
                    #linea = striparray(linea,1)   
                    """
                    pos = findvalinarray(linea,'\n')
                    linea = replaceinarray(linea, pos, '\n')
                    """
                    linea[len(linea)-1] += '\n'
                      
                    filenew.write('\t'.join(linea))
                    i+=1                                  
                    break 
                arrayb = fileb.readline().split('\t') 
    print str(i)+' lines written to '+str(newfile)            
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

def calcmean(line,pos):
    sum = 0    
    failcount=0
    for i in pos:
        try:
            #sum+=float(line[i].replace(',','.'))
            sum+=float(line[i])
        except ValueError:
            failcount+=1

    return sum/(len(pos)-failcount)

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

def shortentitle():
    """
    file = open('CD71GenexpressionsrohdatenNEUGENERIERT.TXT','r')
    filenew = open('CD71GenexpressionsrohdatenNEUGENERIERT_2.TXT','w')
    temp = file.readline()
    temp = temp.split('\t')
    print temp
    for i in range(len(temp)):
        temp[i] = temp[i].replace(' test_(HuEx-1_0-st-v2).rma-exon-all-dabg-Signal','')
        temp[i] = temp[i].replace(' test_(HuEx-1_0-st-v2).rma-exon-all-dabg-Detection ','')
        temp[i] = temp[i].replace('_(HuEx-1_0-st-v2).rma-exon-all-dabg-Signal','')        
        temp[i] = temp[i].replace('_(HuEx-1_0-st-v2).rma-exon-all-dabg-Detection ','')
        print temp[i]    
    temp = '\t'.join(temp)
    print temp
    filenew.write(temp)
    for line in file:
        filenew.write(line)
    """  
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
#getgenpos('CD71GenexpressionsAnnotationBIG.txt',['Gene Symbol','Chromosome','Start']),getgenpos('refGene.txt',['name2','chrom','txStart','txEnd']
def fillemptygenes(arraya,arrayb,xpos,ypos):
    sum = 0.0
    n=0
    for linea in arraya:
        
        if(linea[xpos[0]] == '---'):
            c=0
            for lineb in arrayb:
                if linea[xpos[1]] == lineb[ypos[1]] and ( linea[xpos[2]]  >= lineb[ypos[2]] and linea[xpos[2]]  <= lineb[ypos[3]]):
                  
                    c+=1 
            sum+=c  
            
            n+=1
            if not(n%1000):
                print  'n '+ str(n)+' durchschnitt '+str(sum/n)
                    
    print 'summe : '+str(sum/n)

    return
    

def runit():
    
    meth = readfile('HumanMethylation450_15017482_v.1.2.csv',',',50000)
    genexpr = readfile('CD71Genexpressionsannotation_noemptygenesID.tsv','\t',50000)   
    meth = stripper(meth,getgenpos(meth,'UCSC_RefGene_Name'))  
    meth = expand_and_delete(meth,getgenpos(meth,'UCSC_RefGene_Name'),';')  
    
    #genexpr = readfile('CD71GenexpressionsAnnotationBIG.txt','\t',0)   
    #refGene = readfile('refGene.txt', '\t',0)
    #meth = fillemptygenes(genexpr,refGene,getgenpos(genexpr,['Gene Symbol','Chromosome','Start']),getgenpos(refGene,['name2','chrom','txStart','txEnd']))  
    
    kombi = methy_genexpr_kombinations(meth,getgenpos(meth,'UCSC_RefGene_Name'),genexpr,getgenpos(genexpr,'Gene Symbol'),['UCSC_RefGene_Name','CHR','IlmnID','MAPINFO'],['Probe Set ID','Start'])    
    del meth
    del genexpr
    kombi = sortlist(kombi,getgenpos(kombi, 'IlmnID'),1)
    writelist(kombi,'kombinationen.tsv',0, '\t', 1)
    del kombi
    
    
    #readfile methraw
    #methraw = replacvalueinearray(methraw,',','.')
    #methraw = sortlist(methraw,getgenpos(methraw, 'TargetID'),1)
    #write methraw
    #methraw = readfile('CD71Methyrohdaten_TargetIDsorted.txt', '\t', 0)
    #kombi = appendcolbyarray(kombi,methraw, getgenpos(kombi, 'IlmnID'), getgenpos(methraw, 'TargetID'), getgenpos(kombi, 'UCSC_RefGene_Name'))    
    appendcolbyfile('kombinationen.tsv','CD71Methyrohdaten_TargetIDsorted.txt', getgenpos('kombinationen.tsv', 'IlmnID','\t'), getgenpos('CD71Methyrohdaten_TargetIDsorted.txt', 'TargetID','\t'), getgenpos('kombinationen.tsv', 'UCSC_RefGene_Name','\t'),'kombinationenMethroh.tsv')
    temp = readfile('kombinationenMethroh.tsv', '\t', 150)
    printlist(temp, 0)
    # hier wird die sortierte benutzt evlt vorher sotieren
   
    #genexprraw = readfile('CD71GenexpressionsrohdatenNEUGENERIERT_sortiert_new.TXT', '\t', 0)
    #genexprraw = sortlist(genexprraw,getgenpos(genexprraw, 'Probe Set ID'),1)  
    
    kombi = readfile('kombinationenMethroh.tsv','\t', 0)
    kombi = sortlist(kombi,getgenpos(kombi, 'Probe Set ID'),1)  
    writelist(kombi,'kombinationenMethroh.tsv',0, '\t',0)
    del kombi
    appendcolbyfile('kombinationenMethroh.tsv', 'CD71Genexpr_ProbeSetIDsorted.tsv', getgenpos('kombinationenMethroh.tsv', 'Probe Set ID','\t'), getgenpos('CD71Genexpr_ProbeSetIDsorted.tsv', 'Probe Set ID','\t'), getgenpos('kombinationenMethroh.tsv', 'UCSC_RefGene_Name','\t'),'kombinationenMethrohGenexproh.tsv',1)        
    temp = readfile('kombinationenMethrohGenexproh.tsv', '\t', 150)
    printlist(temp, 0)
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


#def test1(a):
    

#def test2(a):

if __name__ == '__main__':  
    taketime(runit, 1)
    
    
    