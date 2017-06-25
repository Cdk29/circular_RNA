import re
import os

import requests

#this file provide some extra functionnality for control quality for the output report of the circular RNA design process
#here this function allow to trim human microRNAs from the microRNAs of other species, in the file mature.fa, download on microRNA.org
#this microRNAs are then test for alignement on the circular RNA output by the design process 


def mature_reader(file_name):
    #this function is about to trim human microRNAs from microRNAs of other species in the file provide by microRNA.org
    #the name of the file is mature.fa
    #the running time of this function is about few seconds, there is no need to carry an human microRNAs files
  
    ligne=""
    interest=False
    human_microRNAs=[]
    
    mature_miRNAs=open(file_name, "r")
    
    for micro_RNA in mature_miRNAs :
        ligne=micro_RNA
        if interest==True:
            human_microRNAs.append(ligne)
            interest=False
        
        if re.search("Homo", ligne):

            interest=True
            human_microRNAs.append(ligne)
            
    output=open("human_microRNAs.fasta", "w")        
    
    for fasta_line in human_microRNAs :
        output.write(fasta_line) 
        
    output.close()
    mature_miRNAs.close()      

    return 


def miranda_score_reader(output_miranda):
    #the goal of this function is grep the alignements of microRNAs on sequence (binding site or whole circRNA sequence)
    #the function return a list of dictionnary, where dictionnary are on the following form :
    # {'Hit': "   Query:    3' ugcUGCGA-GGGACGACGGUCCa 5'\n                   ::|:| |:|||:||||||| \n   Ref:      5' agaGTGTTGCTCTGTTGCCAGGc 3'\n",
    #'Mir_name': 'hsa-miR-8073',
    #'Score': '166.000000'}]
    
    
    human_microRNAs=open(output_miranda, 'r')
    hits = []
    hit = []
    mir_names = []
    interest = False
    hit_flag=False
    
    
    
    for line in human_microRNAs.readlines():
        if line.startswith("=-"):
            hit = []
            hit_flag=True 
        if line.startswith(" ") and hit_flag==True :
            interest = True
        if interest and line.startswith(" ") and not line.startswith(">hsa"):
            hit.append(line)
        if interest and line.startswith(">hsa"):
            mir_name_holder=line.split('\t')
            interest = False
            mir_name_holder[0]=mir_name_holder[0].replace(">", "")
            mir_names.append(mir_name_holder[0])
            hits.append(hit)
            hit_flag=False
    
    human_microRNAs.close()  
    
    list_of_dictionnary=[]
    
    for hit, mir_name in zip(hits, mir_names):
        
        mir_name_holder=""
        mir_name_holder= str(mir_name)
        
        score_holder=hit[0].split(" ")
        score_holder=score_holder[4]
        
        mir_name={}
        mir_name["Mir_name"]=mir_name_holder
        mir_name["Score"]=score_holder
        mir_name["Hit"]= hit[1] + hit[2] + hit[3]
        list_of_dictionnary.append(mir_name)


    
    
    return list_of_dictionnary






def selecting_best_hit_score(list_of_dictionnary, number_of_best_hits):

    best_hits=[]   #a list of dictionnary of all the better miRanda hits, meaning the alignements with the best score
    
    
    for hits in xrange(0, number_of_best_hits):
        best_score=0
        position_best_mir=0
        for mir_name, x in zip(list_of_dictionnary, xrange(0, len(list_of_dictionnary))):
        
            if mir_name["Score"] >= best_score :
                best_score = mir_name["Score"]
                position_best_mir=x
        
        
        new_list_of_dictionnary=[]
        

        for x in xrange(0, len(list_of_dictionnary)):
            if x==position_best_mir :
                
                best_hits.append(list_of_dictionnary[x])
                continue
                
            new_list_of_dictionnary.append(list_of_dictionnary[x]) 
            
        list_of_dictionnary=new_list_of_dictionnary
    
        
    return best_hits
    





#x="miranda human_microRNAs.fasta hsa_circ_0007874 -noenergy -strict -sc 160 > output_miranda"

#os.system(x) 



#list_of_dictionnary=miranda_score_reader("output_miranda")
#list_of_dictionnary=selecting_best_hit_score(list_of_dictionnary, 3)

















