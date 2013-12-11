from numpy import *
import matplotlib.pyplot as plt
import math

def taketime(f,c,*args):    
    import time 
    start= time.clock()
    for i in range(c):
        f(*args)
    end= time.clock()
    time= (end-start)/c    
    print '\ndurchlaufe: '+str(c) +'\nzeit pro durchlauf: '+str(time) +'\nzeit gesamt: '+str(c*time)
    return time  


def efunktion():
    x = []
    y = []    
    for i in range(25):
        x.append(i)
        y.append(10*(math.exp( -0.3*i ))+2)
    fig = plt.figure(1, figsize=(5,3))
    plt.plot(x,y,color='black',label=r'$10*e^{-0.3*x}+2$')
    #plt.title(r'$\sigma_i=15$')
    #plt.title(r'$10*e^{-0.3*x}+2$')
    plt.legend()
    plt.ylabel('y')
    plt.xlabel('x')
    plt.xlim(0,12)
    plt.ylim(0,12)
    plt.grid(True)
    plt.subplots_adjust(bottom=0.14)
    plt.show()

def genSucheLinear(gen):
    file = open('test/CD71GenexprAnnoNRoh_GeneStartSorted.tsv','r')
    for line in file:
        if line.split('\t')[1]==gen:
            print file.tell()
            return 'found'
    
def genSucheMitIndex(gen):
    f = open('indexfiles/CD71GenexprAnnoNRoh_GeneStartSorted_index.tsv','r')
    off = 0
    f.readline()
    for l in f:
        line = l.strip().split('\t')
        if gen[:2] == line[0]:
            f.close()
            off= long(line[1])
            break
    f.close()
    f = open('test/CD71GenexprAnnoNRoh_GeneStartSorted.tsv','r')
    f.seek(off)
    print off
            
    for line in f:
        line = line.split('\t')
        if line[1] == gen:
            return 'found'
    

#ANLN,LINS,UGCG
taketime(genSucheLinear, 1,'ANTXRL')
taketime(genSucheLinear, 1,'LOC642366')
taketime(genSucheLinear, 1,'TRUB2')
taketime(genSucheMitIndex, 1,'ANTXRL')
taketime(genSucheMitIndex, 1,'LOC642366')
taketime(genSucheMitIndex, 1,'TRUB2')
