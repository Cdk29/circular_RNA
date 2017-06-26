import requests
import tempfile
import zipfile



def RNAfold_latex_output(best_circ, best_score, best_seq):
    
    name_best_cir=str(best_circ)
    name_best_cir=name_best_cir.replace("_", "\_")
    name_best_cir=name_best_cir.replace(">", "$>$")
    best_seq=best_seq.replace("\n", "")
    
    string=""
    string+="\\section{Calculation of the minimum free energy of the circular RNAs}"+ "\n" 
    string+="In this section is report the calculation of the minimum free energy of all the circular RNAs outputed by the design script" +  "\\newline "
    string+="The calculation have been made using RNAfold of the package ViennaRNA 2.3.5. \cite{vienna}" "\\newline "
    string+="Among all the circular RNAs present in the file named \"circular\_design\_output\", the one with the lower minimum free energy is " + name_best_cir + "\\newline "
    string+="The minimum free energy for this design is " + str(best_score) + "kcal/mol" +"\\newline "
    string+="Which sequence is :" + "\\newline " + "\\newline "+ "5'  " + "\\seqsplit{%" + "\n"
    string+= str(best_seq) + "}" + "  3'" + "\n" 

    tex_file = "design_report.tex"
    
    
    fichier=open(tex_file, "a")
    fichier.write(string)
    fichier.close()
    
    return


def latex_output_miranda_header():
    
    string=""
    string+="\\section{Alignements using Miranda}" 
    string+="This section is the report of the best alignements among all the matures human micro-RNAs against the bindings site of micro-RNA on the sponge, and against the whole sequence of the circular RNA produced by the executable." + "\\newline "
    string+="The purpose of this is to check if wether or not the best alignemenents are perform by the micro-RNAs given as arguments to the circular RNA design"
    string+=", and to ensure that a binding site for another micro-RNAs has not been create by mistake somewhere on the circular RNAs sequence during the design process." + "\\newline "
    string+="The list of all the mature micro-RNAs have been download at mirBase.org on May 2017."+ "\\newline "
    string+="The alignements are done using miRanda, with the following options : '-noenergy -strict -sc 150', meaning :" + "\\newline "
    string+="-noenergy : Turn off thermodynamic calculations from RNAlib. If this is used, only the alignment score threshold will be used." + "\\newline "
    string+="-strict : Require strict alignment in the seed region (offset positions 2-8). This option prevents the detection of target sites which contain gaps or non-cannonical base pairing in this region." + "\\newline "
    string+="-sc score : Set the alignment score threshold to score. Only alignments with scores $>=$ score will be used for further analysis." + "\\newline "
    
    
    tex_file = "design_report.tex"
    
    
    fichier=open(tex_file, "a")
    fichier.write(string)
    fichier.close()
    
    return


def latex_output_miranda(list_of_dictionnary, mir):
    
    string=""
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
                holder= "\subsection{" + mir_name + " : " + str(dictionnary[keys])+", cluster for : "  + str(mir) + "}"
                string+= holder
                    
            
            if keys == "Score" :
                    
                    
                holder = "\n" + str(keys) + " : of the alignement :" + str(dictionnary[keys])  + "\n"
                string+=holder
        
    
    
    
    tex_file = "design_report.tex"
    
    
    fichier=open(tex_file, "a")
    fichier.write(string)
    fichier.close()
    
 

    return

def latex_output_cluster(list_of_cluster, list_of_mir, set_of_mir_sequence):
    
    string="""\\newline
   """


    
    for x, mir, seq in zip(xrange(0, len(list_of_cluster)), list_of_mir, set_of_mir_sequence):
     
    
        string+="\\subsection{"+ str(mir)+ "}\n"
        string+="Sequence microRNA : 5' " + str(seq) + " 3' \\newline \\newline"
        string+="3'  " + "\\seqsplit{%" + "\n"
        string+= str(list_of_cluster[x]) + "}" + "  5'" + "\n"


        
    tex_file = "design_report.tex"
    
    fichier=open(tex_file, "a")
    fichier.write(string)
    fichier.close()
    
 

    return

def latex_init():
    
    filin = open("template.txt")
    lines=filin.readlines()
    tex_file = "design_report.tex"  
    fichier=open(tex_file, "w")
    
    for line in lines :

        fichier.write(line)
    fichier.close()

    filin.close()
    
    
    return 


def latex_elong():

    filin = open("template2.txt")
    lines=filin.readlines()
    tex_file = "design_report.tex"  
    fichier=open(tex_file, "a")
    
    for line in lines :

        fichier.write(line)
    fichier.close()

    filin.close()
    
    
    return 


def latex_report_pdf():
    
    #this function use the web service provide by Martin Scharm, called Texpile, available here github.com/binfalse/TEXPILE
    #here the report in latex is rewritte in pdf with the return of texpile

    string=""""	\\bibliography{bibliography}
    \\end{document} 
    """
        
    
    tex_file = "design_report.tex"
    bib_file = "bibliography.bib"
    
    fichier=open(tex_file, "a")
    fichier.write(string)
    fichier.close()
    

    tmp = tempfile.NamedTemporaryFile(delete=False)    
    z = zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED)
    z.write(tex_file)
    z.write(bib_file)
    z.close()
    
    print (tmp)


    # send the project to TEXPILE
    r = requests.post('http://texpile.bio.informatik.uni-rostock.de', files={'project': open(tmp.name, 'rb')}, data={"filename": tex_file})

    
    
    
    #r = requests.post('http://texpile.bio.informatik.uni-rostock.de', files={'project': open(tex_file, 'rb')})
    with open(tex_file + ".pdf", 'wb') as result:
        result.write(r.content)

   
    
    return 

