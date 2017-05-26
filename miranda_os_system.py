import re
import os


#miranda for the circular RNA report outpout





def mature_reader(file_name):
    #this function is about to trim human microRNAs from microRNAs of other species in the file provide by microRNA.org
    #the name of the file is mature.fa
    #the running time of this function is about few seconds, there is no need to carry an human microRNAs files
  
    ligne=""
    interest=False
    human_microRNAs=[]
    
    mature_miRNAs=open("mature.fa", "r")
    
    for micro_RNA in mature_miRNAs :
        ligne=micro_RNA
        if interest==True:
            human_microRNAs.append(ligne)
            interest=False
        
        if re.search("Homo", ligne):
            print "ok"
            interest=True
            human_microRNAs.append(ligne)
            
    output=open("human_microRNAs.fasta", "w")        
    
    for fasta_line in human_microRNAs :
        output.write(fasta_line) 
        
    output.close()
    mature_miRNAs.close()      

    return 




x="miranda mature.fa hsa_circ_0007874 -noenergy -strict -sc 160 > output_miranda"



os.system(x)   #terminal is "froze" the time to execute the query, thanksfully 


filin=open("mature.fa", "r")

read=filin.readlines()

for ligne in read :
    print "line", ligne
    































