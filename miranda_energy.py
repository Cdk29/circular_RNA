
# miranda human_microRNAs.fasta hsa_circ_0007874 -strict > output_miranda_energy 

def miranda_energy_reader(output_miranda_energy):
    #the goal of this function is grep the alignements of microRNAs on sequence from the circular RNA, based on their energy
    #the function return a list of dictionnary, where dictionnary are on the following form :

    
    
    human_microRNAs=open(output_miranda_energy, 'r')
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
        
        energy_holder=hit[4].split(" ")
        energy_holder=energy_holder[5]
        energy_holder=float(energy_holder)

        mir_name={}
        mir_name["Mir_name"]=mir_name_holder
        mir_name["Energy"]=energy_holder
        mir_name["Hit"]= hit[1] + hit[2] + hit[3]
        list_of_dictionnary.append(mir_name)

    
    
    
    return list_of_dictionnary



def selecting_best_hit_energy(list_of_dictionnary, number_of_best_hits):

    best_hits=[]   #a list of dictionnary of all the better miRanda hits, meaning the alignements with the best score
    
    
    for hits in xrange(0, number_of_best_hits):
        lowest_energy=0
        position_best_mir=0
        for mir_name, x in zip(list_of_dictionnary, xrange(0, len(list_of_dictionnary))):

            if mir_name["Energy"] < lowest_energy :
                lowest_energy = mir_name["Energy"]
                position_best_mir=x
        
        
        new_list_of_dictionnary=[]
        

        for x in xrange(0, len(list_of_dictionnary)):
            if x==position_best_mir :
                
                best_hits.append(list_of_dictionnary[x])
                continue
                
            new_list_of_dictionnary.append(list_of_dictionnary[x]) 
            
        list_of_dictionnary=new_list_of_dictionnary
    
        
    return best_hits






