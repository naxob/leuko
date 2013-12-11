fmyth = open('CD71MethylierungsAnnotation_unique_noemptygenesID.tsv','r')
fgen = open('CD71Genexpressionsannotation_noemptygenesID.tsv','r')

fmyth.readline()
fgen.readline()

myth = fmyth.readlines()
gen = fgen.readlines()

dmyth = {}
dgen = {}
for x in myth:
    temp = x.split('\t')
    dmyth.setdefault(temp[4])

for y in gen:
    temp =y.split('\t')
    dgen.setdefault(temp[2])

print 'anzahl gene in myth : '+str(len(dmyth))
print 'anzahl gene in genexpr : '+str(len(dgen))

fgen.close()
fmyth.close()