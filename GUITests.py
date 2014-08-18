'''
Created on 30.07.2014

@author: naxobIdeaPad
'''
import numpy as np
import os

if __name__ == '__main__':
    a = []
    
    for file in os.listdir("C:/Users/dvrb2/Downloads/WIN98Virtualisierung/ImaGene/v4/ImaGene40/JavaSources/Jarextract/all/"):
        if file.endswith(".java"):
            a.append(file)
    print a
      
    fw = open('C:/Users/dvrb2/Downloads/WIN98Virtualisierung/ImaGene/v4/ImaGene40/JavaSources/Jarextract/all/erg/search.txt','w')     
    
    for path in a:
        f = open('C:/Users/dvrb2/Downloads/WIN98Virtualisierung/ImaGene/v4/ImaGene40/JavaSources/Jarextract/all/'+path,'r')
        for line in f:
            if line.find('config') != -1:
                print path+' '+line
                fw.write(path+'\n')
            
        f.close()
    fw.close()
    