fGenexpAnno = open('CD71Genexpressionsannotation_noemptygenes.tsv','r')
fMythAnno = open('CD71MethylierungsAnnotation_unique_noemptygenes.tsv','r')

print '\n Genexpression -------------- \n'

out = 0
for y in fGenexpAnno :
    gen_row = y.split('\t') 
    print gen_row[1]
            
    inner = 0        
    for x in fMythAnno :
        myth_row = x.split('\t') 
        # print 'Mythelierung : '+myth_row[4]
        if myth_row[4]==gen_row[1] : 
            print '\n -----  MATCH  ----- \n'
        inner+=1
        if inner>1000 : break 
    print out
    out+=1    
    if out>20 : break
    
print '----  done'
fGenexpAnno.close
fMythAnno.close