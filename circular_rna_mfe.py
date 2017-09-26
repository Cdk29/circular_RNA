import os

#this set of function are used to calculate the minimum free energy structure of several circular RNA outputed by the circularRNA
#design script, and  select the circular RNA with the less Minimum Free Energy structure,


#considering the use of viennaRNA package, and more specifically, the use of the RNAfold executable, the output if the circularRNA design
#script will be a text file containing several sequence separate by a line starting with a ">" caractec and the name of the circular RNA

def ViennaRNA_RNAfold(circular_design_output):
    #this function take the output of the circularRNA design amd
    
    print "RNAfold is running on the output of the circular RNA design script"
    
    command="RNAfold --circ " + circular_design_output + " > outputvienna"
    os.system(command)


    return


def RNAfold_grep():


    
    fichier=open("outputvienna", "r")

    best_score=100
    
    for line in fichier.readlines():
        
        line=line.replace("( -", "(-")
        line=line.replace("  0.00", "0.00")
        split=line.split()

        if line.startswith(">"):
            circ_name=line
            
        if len(split)==1:
            sequence=line

        if len(split)==2 :

            score = split[1]
            score = score.replace("(", "")
            score = score.replace(")", "")
            score = float(score)

            if score<best_score:

                best_score=score
                best_seq=sequence
                best_circ=circ_name
        
    
    fichier.close()

    return best_circ, best_score, best_seq




#circular_design_output="test_file_for_vienna"


#ViennaRNA_RNAfold(circular_design_output)


#best_circ, best_score, best_seq = RNAfold_grep()

































