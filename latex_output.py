import re
import os
import requests


    

def latex_output_miranda(list_of_dictionnary):
    
    string=""
    string+="\\section{Miranda Alignements}" 
    string+="This section is the report of the best alignements among all the matures human micro-RNAs against the bindings site of micro-RNA on the sponge, and against the whole sequence of the circular RNA produced by the executable." + "\\newline "
    string+="The purpose of this is to check if wether or not the best alignemenents are perform by the micro-RNAs given as arguments to the circular RNA design"
    string+=", and to ensure that a binding site for another micro-RNAs has not been create by mistake somewhere on the circular RNAs sequence during the design process." + "\\newline "
    string+="The list of all the mature micro-RNAs have been download at mirBase.org on May 2017."+ "\\newline "
    string+="The alignements are done using miRanda, with the following options : '-noenergy -strict -sc 160', meaning :" + "\\newline "
    string+="-noenergy : Turn off thermodynamic calculations from RNAlib. If this is used, only the alignment score threshold will be used." + "\\newline "
    string+="-strict : Require strict alignment in the seed region (offset positions 2-8). This option prevents the detection of target sites which contain gaps or non-cannonical base pairing in this region." + "\\newline "
    string+="-sc score : Set the alignment score threshold to score. Only alignments with scores $>=$ score will be used for further analysis." + "\\newline "
    
    for dictionnary in list_of_dictionnary:
        
        for keys in dictionnary:
            if keys == "Hit" :
                    
                holder=dictionnary[keys]
                    
                holder=holder.replace("Query:", "miRNA:")
                holder=holder.replace("Ref:    ", "CircRNA:")
                
                holder= "\\begin{verbatim}" + "\n" + holder
                holder=holder +"\n" + "\end{verbatim}"
                string+=holder
    
            if keys == "Mir_name" :
                
                mir_name=str(keys)
                mir_name=mir_name.replace("_", "\_")
                holder= "\subsection{" + mir_name + " : " + str(dictionnary[keys])+"}"
                string+= holder
                    
            
            if keys == "Score" :
                    
                    
                holder = "\n" + str(keys) + " : of the alignement :" + str(dictionnary[keys])  + "\n"
                string+=holder
        
    
    
    
    tex_file = "design_report.tex"
    
    fichier=open(tex_file, "a")
    fichier.write(string)
    fichier.close()
    
 

    return



def latex_report_pdf():
    
    #this function use the web service provide by Martin Scharm, called Texpile, available here github.com/binfalse/TEXPILE
    #here the report in latex is rewritte in pdf with the return of texpile

    string="\\end{document}"
        
    
    tex_file = "design_report.tex"
    
    fichier=open(tex_file, "a")
    fichier.write(string)
    fichier.close()
    
    
    r = requests.post('http://texpile.bio.informatik.uni-rostock.de', files={'project': open(tex_file, 'rb')})
    with open(tex_file + ".pdf", 'wb') as result:
        result.write(r.content)

   
    
    return 

