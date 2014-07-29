import File
import Plots
'''
Created on 18.07.2014

@author: Axel Wilbertz
'''

if __name__ == '__main__':
    
    file1 = File.FileObject("Cor_CD71MethExpCore_R_TOP10k.tsv",0,'\t')
    file1.getGeneValues("ARV1",False)
    file1.getTitles()
    mapper = File.Mapper('maps\KorrelationMapping.csv',':');
    mapper.extractGroupMapping()
    print 'titles: '+str(file1.titles)
    print 'values: '+str(file1.values)
    print 'mapper: '+str(mapper.mapperdict)
    mapper.findPositions(file1.titles,mapper.mapperdict)
    print 'pos: '+str(mapper.positions)
    corplot = Plots.CorrelationPlot(file1.titles,file1.values,file1.genecolumn,['g','r','y','c','b'],mapper.positions);
    corplot.plotValues()

    print '------------------ splicing ---------------'
    file2 = File.FileObject("EXON-EXTENDEDshort.tsv",1,'\t')
    file2.getGeneValues("TP53",True)
    file2.getTitles()
    mapper2 = File.Mapper('maps\SplicingTranscriptExpressionMapping.csv',':');
    mapper2.extractGroupMapping()
    mapper2.findPositions(file2.titles,mapper2.mapperdict)
  
    file3 = File.RefGene('refgene/refGene.tsv','indexfiles/refGene_index.tsv',file2.genecolumn,'\t',1,9,10,7,8,file2.titles,file2.values)
    file3.extractRefGeneIndexFileEntries()
    file3.extractRefGeneValues()
    file3.calcIntronExon()    
      
    splicingplot= Plots.SplicingPlot(file2.titles,file2.values,file2.genecolumn,['g','r','y','c','b'],mapper2.positions,file3.intron,file3.exon)
    splicingplot.setMDSSubtypes()
    splicingplot.calcMeanANDStd()
    splicingplot.plotValues()

    
    
    

