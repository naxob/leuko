"""print 'working'
fmeth = open('HumanMethylation450_15017482_v.1.2.csv','r')
methtitle = fmeth.readline().split(',')
methvalues = fmeth.readline().split(',')

i=0
for x in methtitle:
    print str(i)+' '+x+'\t\t'+methvalues[i]  
    i+=1


    
"""
print 'working'
fmeth = open('HumanMethylation450_15017482_v.1.2.csv','r')
meth = fmeth.readlines()
metharray= []
for x in meth:
    metharray.append(x.split(','))

i=0
for x in metharray:
    if x[0] != x[1]:
        print str(i)+' '+x[0]+' not equal '+x[1]
        if x[0].strip() == x[1].strip:
            print 'but equal after strip'
    i+=1

fmeth.close()
print 'done'