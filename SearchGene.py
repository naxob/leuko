print 'started'
# fnew = open('CD71Genexpressionsannotation_noemptygenesID.tsv','r')
fnew = open('CD71GenexpressionsAnnotationBIG.txt','r')

genfile = fnew.readlines()
fnew = open('CD71MethylierungsAnnotation_unique_noemptygenesID.tsv','r')
mythfile = fnew.readlines()
k=0
for x in genfile:
    temp = x.split('\t')
    if temp[1] == 'PITX1':
            print x
            k+=1
print '---------- anzahl '+str(k)    
l=0       
for x in mythfile:
    temp = x.split('\t')
    if temp[4] == 'PITX1':
            print x
            l+=1
print '---------- anzahl '+str(l)             
print 'moegliche Kombinationen '+str(l*k)
print 'done'
 
fnew.close()
"""
print 'start'
fmyth = open('CD71MethylierungsAnnotation_unique_noemptygenesID.tsv','r')
fmyth.readline()
myth = fmyth.readlines()

fgen = open('CD71Genexpressionsannotation_noemptygenesID.tsv','r')
fgen.readline()
gen = fgen.readlines()

m={}
for x in myth:
    temp = x.split('\t')
    m.setdefault(temp[4])
print len(m)

g={}
for x in gen:
    temp = x.split('\t')
    g.setdefault(temp[2])
print len(g)

sm=set(m)
sg=set(g)
ss = sm&sg
sd1 = sm-sg
sd2 = sg-sm
print 'laenge myth : '+str(len(sm))
print 'laenge gen : '+str(len(sg))
print 'sd1 : '+str(len(sd1))
print 'sd2 : '+str(len(sd2))

fnew = open('NichtSchnittmengeGenexprUndMyth.txt','w')
fnew.write('not in genexpr\n')
for x in sd1:
    fnew.write(x+'\n')
fnew.write('not in myth\n')
for x in sd2:
    fnew.write(x+'\n')
fnew.close()
"""
print 'done'