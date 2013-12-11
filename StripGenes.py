print 'working...'
fnew = open('CD71Genexpressionsannotation_noemptygenesID.tsv','r')
fnew2 = open('CD71Genexpressionsannotation_noemptygenesID2.tsv','w')
fnew2.write(fnew.readline())
file = fnew.readlines()
genarray = []
for x in file:
    genarray.append(x.split('\t'))  

for x in genarray:
    temp = x[2].strip()
    x[2] = temp
    fnew2.write('\t'.join(x))

print '...done'