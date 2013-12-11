import sets
f = open('CD71MethylierungsAnnotation_tabstop.txt','r')
fnew = open('writetest.csv','w')

i = 0
for line in f:    
    if i == 0 : fnew.write(line)
    else :
        # Trennt die Spalten nach tabstop 
        linesplit = line.split('\t')
        # wenn kein gen eingetragen in neue file schreiben
        if linesplit[4] == '' :
            # fnew.write(line)
            print 'i: ' + str(i) +' write ' + ' copy'
        # mehrere gene eingetragen duplikate entfernen dann in file schreiben
        elif linesplit[4].count(';')>=1 :   
            print 'i: ' + str(i) + '\t value: ' + linesplit[4] + '\t\t\t\t! edit !'    
            # vorne und hinten wurden """ angefuegt danke excel
            temp = linesplit[4].split('"""')
            
            colsplit_duplicates = temp[1].split(';')
            # entferne doppelten gene aus der liste durch wandeln in set
            colsplit_unique = list(sets.Set(colsplit_duplicates))
            print colsplit_unique   
            # wenn nur noch ein gen vorhanden ist schreiben in file
            if len(colsplit_unique)==1 :
                linesplit[4]=colsplit_unique[0]
                wtemp = '\t'.join(linesplit)
                fnew.write(wtemp)
            # wenn mehr als ein gen vorhanden ist dann muss zeile kopiert werden und neues gen eingefuegt
            else :
                for x in colsplit_unique :
                    linesplit[4]=x
                    wtemp = '\t'.join(linesplit)
                    fnew.write(wtemp)
            
        # wenn nur ein gen eingetragen ist in file schreiben
        else :
            fnew.write(line)
            print 'i: ' + str(i) + '\t value: ' + linesplit[4] +' single gene'        
            
    i+=1
    if i > 1000000 : break    
    
f.close()
fnew.close()