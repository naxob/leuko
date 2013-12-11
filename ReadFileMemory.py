print 'working...\n'
fm = open('CD71MethylierungsAnnotation_unique_noemptygenesID.tsv','r')
fnew = open('kombinationen.tsv','w')
fnew.write('Gene\tTargetID\tProbeSetID\tChrom\tMMapinfo\tGStart\n')
# erste Zeile wegschmeissen
fm.readline()
fmm = fm.readlines()
fm.close()
fg = open('CD71Genexpressionsannotation_noemptygenesID.tsv','r')
# erste Zeile wegschmeissen
fg.readline()
fgm = fg.readlines()
fg.close()

print 'Files in Speicher gelesen'


# liest file in array ein

mytharray = []
i=0
for x in fmm:
    mytharray.append(x.split('\t'))    
    if i > 10000000 : break
    i+=1
    
# liest file in array ein
genarray = []
y=0
for x in fgm:
    genarray.append(x.split('\t'))    
    if y > 100000000 : break
    y+=1
print 'Files in Arrays geschrieben und nach Tabstop getrennt'

mythdict = {}
for x in mytharray :       
    mythdict.setdefault(x[4], []).append(x[0])

gendict = {}
for x in genarray :
    gendict.setdefault(x[2], []).append(x[0])

print 'Dictionary gefuellt' 
wfile = open
while(len(mythdict)):
    # print 'laenge myth\t'+str(len(mythdict))+'\tlaenge gen\t'+str(len(gendict))
    mythgen = mythdict.popitem()    
    try :
        gengen = gendict.pop(mythgen[0])
        for x in mythgen[1]:
            for y in gengen:
                #             Gene            TargetID                     ProbeSetID               MChrom                     MMapinfo                    GStart
                fnew.write(mythgen[0] +'\t'+mytharray[int(x)][1]+ '\t'+genarray[int(y)][1]+'\t'+mytharray[int(x)][2]+'\t'+mytharray[int(x)][3]+'\t'+genarray[int(y)][15]+'\n')
                
    except KeyError:
        pass
print 'Kombinationen in neue Datei geschrieben'       

print '...done'

fnew.close()
