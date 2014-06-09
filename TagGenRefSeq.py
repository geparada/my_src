import sys
import csv


RNA = []
Genoma = []
SeqName = []
Resultados=[]


def SJtag(AliSeq, len_tag, GenomicTags=""):
    '''Generates tags of the splice-juntions
    *Input file = USCS BLAT alignments table (22 columns ) joined with fasta converted to tab (3 columns). This file has 25 columns.
    *Usage = python Tag-Fulinfo inputfile tag length
    *Output = seq ID , intron name, strand, tag sequence, tag length '''
       
    
    reader1 = csv.reader(open(AliSeq), dialect='excel-tab' )
         
         # row[0]  matches
         # row[1]  misMatches
         # row[2]  repMaches
         # row[3]  nCount
         # row[4]  qNumInsert    
         # row[5]  qBaseInsert
         # row[6]  tNumInsert
         # row[7]  tBaseInsert
         # row[8]  strand
         # row[9] qName
         # row[10] qSize
         # row[11] qStart
         # row[12] qEnd         
         # row[13] tName
         # row[14] tSize
         # row[15] tStart
         # row[16] tEnd
         # row[17] blockCount
         # row[18] blockSizes
         # row[19] qStarts
         # row[20] tStarts

         # row[21] qName
         # row[22] 1
         # row[23] sequence

    
   
    for row in reader1:
        
        qstarts = map (int, row[19].split(",")[1:-1])                      # obteniendo lista int de las distintas variables
        tstarts = map(int, row[20].split(",")[0:-1])
        blocksizes = map(int, row[18].split(",")[0:-1])
       
        n = len_tag/2                                                      #preparando argumento para su uso
      
        
        for t, b, t2, b2, q in zip(tstarts, blocksizes, tstarts[1:], blocksizes[1:]  , qstarts):
            
            if b>=0 and b2>=0:                                   #filtro por tamano de exones
                tbsum = t+b                                      # generando primera coordenada del intron al cual se le va a hacer el tag              
                start = q-n                                     #generando intervalos para hacer el slicing del la secuencia                      
                end = q+n
    
                    
                if row[8] == '+' :                          #Para los que aliniean en la hebra +
    
                    
                    if start<0:                             #Precausiones generar buenos tag del primer y ultimo tag
                        start = 0  
                    if end>len(row[23]):
                        end=len((row[23]))
                    tag = row[23][start:end]
    
                    
                if row[8] == '-' :                                #Para los que alinian en la hebra - es todo al inverso
    
                    if end>len(row[23]):
                        end=len(row[23])                                       
                    tag = row[23][-end:-start]
                    if start<=0: 
                        tag = row[23][-end:]
         

                if GenomicTags!="tags":                      
                    name_seq_tag = row[9], row[13] + str(tbsum) + row[8] + str(t2), tag
                    SeqName.append(name_seq_tag)
                if GenomicTags=="tags":
                    #print row[13] + str(tbsum) + row[8] + str(t2) + row[9], row[13] + str(tbsum) + row[8] + str(t2) + tag
		    #print row[13] + str(tbsum) + row[8] + str(t2), tag, len(tag)
		    #print ">" + row[9] + row[13] + str(tbsum) + row[8] + str(t2) + "\n" + tag                
		    print row[13] + str(tbsum) + row[8] + str(t2) + "\t" + row[9]

def switch (file2=""):
    if file2=="tags":
        pass
    else:
        CompareTags(sys.argv[3], sys.argv[4])

           


def CompareTags(GalaxyGenomico, output="F"):

    reader2 = csv.reader(open(GalaxyGenomico), dialect='excel-tab' )
        # row [0] intronName
        # row [9] seq
           
    for intro in reader2:        
        B = intro[0], intro[9]
        Genoma.append(B)

    set_Genoma = set(Genoma)
    result = [x for x in SeqName if (x[1], x[2]) in set_Genoma]
    for c, d, e in result:
        if output=="F":
            print c, d, e
        if output=="D":
            filtrados = c, d, e
            Resultados.append(filtrados)
        
    if output=="D":
        sacados = set(SeqName)-set(Resultados)
        for c, d, e in sacados:
            print c, d, e
               
    


    
    
if __name__ == '__main__':                                      #Permite tomar argumentos mediante el terminal
    SJtag(sys.argv[1], int(sys.argv[2]), sys.argv[3])
    switch(sys.argv[3])

#SJtag("BLAToincDNAmm9",  30, "BLATIntronescDNAmm9")
#switch("BLATIntronescDNAmm9")
#CompareTags("BLATIntronescDNAmm9")
