import RNA
import os


#this set of function are used to calculate the minimum free energy structure of several circular RNA outputed by the circularRNA
#design script, and  select the circular RNA with the less Minimum Free Energy structure,


#considering the use of viennaRNA package, and more specifically, the use of the RNAfold executable, the output if the circularRNA design
#script will be a text file containing several sequence separate by a line starting with a ">" caractec and the name of the circular RNA

def ViennaRNA_RNAfold(circular_design_output):
    #this function take the output of the circularRNA design amd
    
    
    print "RNAfold is running on the output of the circular RNA design script"
    command="RNAfold<"+ circular_design_output + " --circ > outputvienna"
    os.system(command)

    best_circ, best_score, best_seq = RNAfold_grep()


    return


def RNAfold_grep():


    fichier=open("outputvienna", "r")
    
    
    best_score=0
    
    for line in fichier.readlines():
        
        split=line.split()
        
        if line.startswith(">"):
            circ_name=line
            continue
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




circular_design_output="test_file_for_vienna"


ViennaRNA_RNAfold(circular_design_output)




def RNAfold_latex_output(best_circ, best_score, best_seq):
    
    name_best_cir=str(best_circ)
    name_best_cir=name_best_cir.replace("_", "\_")
    name_best_cir=name_best_cir.replace(">", "$>$")
    best_seq=best_seq.replace("\n", "")
    
    string=""
    string+="\\section{Calculation of the minimum free energy of the circular RNAs}"+ "\n" 
    string+="In this section is report the calculation of the minimum free energy of all the circular RNAs outputed by the design script" +  "\\newline "
    string+="The calculation have been made using RNAfold of the package ViennaRNA 2.3.5." "\\newline "
    string+="Among all the circular RNAs present in the file named \"circular\_design\_output\", the one with the lower minimum free energy is " + name_best_cir + "\\newline "
    string+="The minimum free energy for this design is " + str(best_score) + "kcal/mol" +"\\newline "
    string+="Which sequence is :" + "\\newline " + "\\newline "+ "5'  " + "\\seqsplit{%" + "\n"
    string+= str(best_seq) + "}" + "  3'" + "\n" 

    tex_file = "design_report.tex"
    
    
    fichier=open(tex_file, "a")
    fichier.write(string)
    fichier.close()
    
    return




best_circ, best_score, best_seq = RNAfold_grep()



RNAfold_latex_output(best_circ, best_score, best_seq)


































