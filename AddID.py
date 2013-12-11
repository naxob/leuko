
fw = open('CD71MethylierungsAnnotation_unique_noemptygenes.tsv','r')
fnew = open('CD71MethylierungsAnnotation_unique_noemptygenes_2.tsv','r')

i = 0
fnew.write('ID'+fw.readline()) 

for x in fw :    
    temp = x.split('\t')
    temp[0]=str(i)
    fnew.write('\t'.join(temp))    
    i+=1
    
fw.close()
print 'done'
fnew.close()
