import json
import os
import sys
import getopt
import shutil
import requests
from StringIO import StringIO
from lxml import etree
from bs4 import BeautifulSoup
import itertools
from circular_rna_mfe import *
from miranda_score import *
from latex_output import *



distance=17             #this is a global count to separate the beginning of a seed to the beginning of another one from 17 nt
list_of_microRNAs=""
priority=""
size=300
list_of_sequence=""
override=[]

def micro_RNA_sequence_request(mirna):
    
    Accession_ID=triplex_RNA_database_request(mirna)
    micro_RNA_sequence=mirbase_request(Accession_ID)
    micro_RNA_sequence=micro_RNA_sequence.lower()
    return micro_RNA_sequence

def triplex_RNA_database_request(mirna):
    # request to get the Accession ID to iterrogate mirbase from the name of the microRNA
    # query the TriplexRNA for a given miRNA name, and parse the returned
    # HTML tree in the response to get the correspondent miRNA's MIMAT ID
    # DTD for parsing the HTML tree
    dtd = "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">\n"

    query = "https://www.sbi.uni-rostock.de/triplexrna/Human/mirna/" + mirna
    response = requests.get(query)
    if response.status_code == 200:
        tree=parsing_HTML_tree(response, dtd)
        target=xpath(tree, mirna)
        # the result is a list of elements
        if len(target) != 0: 
              Accession_ID=get_mimat_id(target)
    return Accession_ID


def mirbase_request(Accession_ID):
    #request of mirbase with the Accession ID
    response = requests.get("http://www.mirbase.org/cgi-bin/get_seq.pl?acc=" + Accession_ID)
    if response.status_code == 200:
        micro_RNA_sequence=get_micro_RNA_sequence(response)
    return micro_RNA_sequence


def get_micro_RNA_sequence(response):
    
    data_from_mirbase=response.content
    data_list=data_from_mirbase.split("\n")
    micro_RNA_sequence=data_list[2]
    
    return micro_RNA_sequence

def parsing_HTML_tree(response, dtd):
    # parse the HTML tree
    data = BeautifulSoup(dtd + response.content, 'lxml')
    parser = etree.XMLParser(remove_blank_text=True, resolve_entities=False)
    tree = etree.parse(StringIO(data), parser)
    
    return tree


def xpath(tree, mirna):
    # get the element whose text contains the given miRNA name
    
    target = tree.xpath("//table/tbody/tr/td/a[text() = \"" + mirna + "\"]")
    
    return target



def get_mimat_id(target):
    # get the URL in the href attribute, and return the MIMAT ID
    mimat = target[0].xpath('./@href')
    mimat_id = mimat[0].split('=')[1]
    
    return mimat_id


def elongation(sequence, micro_RNA_sequence):
    
    sequence = micro_RNA_binding_site(micro_RNA_sequence, sequence)
    sequence = seed_site_distance(micro_RNA_sequence, sequence)
    
    return sequence


def micro_RNA_binding_site(micro_RNA_sequence, sequence):
    # sequence construction for the binding site of a microRNA 
    # the microRNA sequence is write as 5' to 3'
    # the circularRNA is write as 3' to 5'
    
    sequence=A_anchor_1(sequence)   
    sequence=seed_binding_site(sequence, micro_RNA_sequence)    
    sequence=A_anchor_9(sequence)    
    sequence=clivage_protection(sequence, micro_RNA_sequence)    
    sequence=supplementary_site(sequence, micro_RNA_sequence)    
    sequence=end_processing(sequence, micro_RNA_sequence)

    return sequence

def A_anchor_1(sequence):
    # Conserved seed pairing, often flanked by adenosines, indicates that thousands of human genes are microRNA targets.
    # 10.1016/j.cell.2004.12.035

    sequence = sequence + "A"

    return sequence


def seed_binding_site(sequence, micro_RNA_sequence):
    # creation of the complementary site for the micro_RMA seed 
    # MicroRNAs: Target Recognition and Regulatory Functions DOI 10.1016/j.cell.2009.01.002
    # MicroRNA Targeting Specificity in Mammals:Determinants beyond Seed Pairing DOI 10.1016/j.molcel.2007.06.017

    for i in micro_RNA_sequence[1:8]:
        if i=="a":
            sequence = sequence + "U"
        elif i=="c":
            sequence = sequence + "G"
        elif i=="g":
            sequence = sequence + "C"
        elif i=="u":
            sequence = sequence + "A"
    return sequence 


def A_anchor_9(sequence):
    # Conserved seed pairing, often flanked by adenosines, indicates that thousands of human genes are microRNA targets.
    # 10.1016/j.cell.2004.12.035

    sequence = sequence + "A"

    return sequence

def clivage_protection(sequence, micro_RNA_sequence):
    #avoidance of complementarity between the binding site and the microRNA for the microRNA's nucleotide 9 to 12
    for i in micro_RNA_sequence[9:12]:
        if i == "a" :
            sequence=sequence+"A"
        else:
            sequence=sequence+"X"

    return sequence

def supplementary_site(sequence, micro_RNA_sequence):
    for i in micro_RNA_sequence[12:16]:
        if i=="a":
            sequence = sequence + "U"
        elif i=="c":
            sequence = sequence + "G"
        elif i=="g":
            sequence = sequence + "C"
        elif i=="u":
            sequence = sequence + "A"
    return sequence 



def end_processing(sequence, micro_RNA_sequence):
    #from nucleotide 17 until the end of the microRNA, addition of a A or U with avoidance of complementarity on the binding site
    for i in micro_RNA_sequence[16:]:
        if i=="a":
            sequence=sequence+"A"
        else:
            sequence=sequence+"X"
            
    return sequence




def seed_site_distance(micro_RNA_sequence, sequence):
    # Saetrom et al,(2007) Distance constraints between microRNA target sites dictate efficacy and cooperativity. 
    # DOI 10.1093/nar/gkm133
    # here the choose is to put 17 nucleotide between the beginning of the first seed and the second seed, according
    # to the figure 2.B from Saetrom et al, 2007
    # to get a good comprhension of the structure of the separation, you can refer to the pdf output or fig 2.A from the same article
    
    global distance

    seed_site_distance=len(micro_RNA_sequence[1:len(micro_RNA_sequence)])

    if seed_site_distance < distance :
        while seed_site_distance != distance:
            sequence=sequence+"X"
            seed_site_distance=seed_site_distance+1
    return sequence


def poly_U_liker(sequence):
    #micro-RNA binding site with hight local content of  A and U perform the best 
    #MicroRNAs: Target Recognition and Regulatory Functions
    #DOI 10.1016/j.cell.2009.01.002
    #the goal is here to have a binding site with little potential for forming self-interacting secondary structure
    #there is some protein which seems to target AU rich sequence (AREs sequence), however it seems that they require some A 
    # so we choose here to put only U for a U rich design
    # Spasic et al (2012). Genome-Wide Assessment of AU-Rich Elements by the
    #AREScore Algorithm. PLoS Genetics, 8(1), e1002433. http://doi.org/10.1371/journal.pgen.1002433
    
    
    returned_sequence=""
    for i in sequence:
        if i=="X" or i=="W":

            returned_sequence=returned_sequence+"U"

        else:
            returned_sequence=returned_sequence+i

    return returned_sequence


def from_5_prime_to_3_prime_correction(sequence):
    #the circular RNA is writte from 3 prime to 5 prime, since the microRNA is writte 5 prime to 3 prime
    #this rewritte the sequence in the 5' to 3' direction
    returned_sequence=""
    for i in sequence:
        returned_sequence=i+returned_sequence
    
    return returned_sequence



def print_help():
    
    #Prints the help/usage
    
    print('\t This function, ' + os.path.basename(sys.argv[0]) + ', take as arguments:' + '\n')
    print('\t-h  | --help\t\t\tprint this help and exit\n')
    print('\t-l  | --list_of_microRNAs= \t the name of the file with the microRNAs name, as they are inside the triplexRNA database. \n \t\t\t\t\t One microRNA name by line. Please type the name of the file whitout "" \n')
    print('''\t-p  | --priority= \t\t the name of the file with the priority for each microRNAs, 
               \t\t\t\t in the same order than the microRNAs in the files specifying for the microRNAs names \n \t\t\t\t\t One int per line \n''')
    print('\t-s  | --size= \t   \t\t The size in nucleotides for the circular RNA. Default value : 300 " \n')
    print('\t-d  | --distance=\t\t The distance between seed as define in Saetrom et al,(2007) Distance constraints between microRNA target sites dictate efficacy and cooperativity. Default value : 17 nt \n')
    print('''\t-q  | --list_of_sequence= \t The name of a file with a list of sequence for binding site for microRNAs, in case of problem with the quality control with a previous running \n
         \t\t\t\t This is an override of the normal process to create binding site, in case the quality control reveal that another microRNA is susceptible to bind to the binding site.
         \t\t\t\t This allow to rerun the script with the new binding site to submit it to the quality control.
         \t\t\t\t The sequences has to be passed, one per line, in the same order than the microRNAs that they target, in the file listing the microRNAs.
         \t\t\t\t To override the process for only one microRNA, leave the other line empties.
         \t\t\t\t This sequences have to been write from 3' to 5'.
          ''')
    print('\tPlease note the "=" after the longs options. The short arguments don\'t take "=", juste a space after -x')
    print('\n' + '\t' + 'The script proceed by extended every cluster of seed by a number of iteration based on the priority number')


    return




def get_cli(argv):
    """
    Parses the command line, returning the mandatory parameters.
    """

    global list_of_microRNAs
    global priority  
    global size
    global list_of_sequence
    global distance

    try:
        opts, args = getopt.getopt(argv, "hl:p:s:q:d:", ['help', 'list_of_microRNAs=', 'priority=', 'size=', 'list_of_sequence=', 'distance=' ])
    except getopt.GetoptError as err:
        print_help()
        print str(err)
        sys.exit(1)
    
    # parse the provided cli options
    
    for opt, arg in opts:
        print opt, arg
        if opt in ('-h', '--help'):
            print_help()
            sys.exit(0)
        elif opt in ( '-l', '--list_of_microRNAs'):
            list_of_microRNAs = arg
        elif opt in ('-p', '--priority'):
            priority = arg
        elif opt in ('-s', '--size'):
            size = int(arg)
        elif opt in ('-q', '--list_of_sequence='):
            list_of_sequence = arg
        elif opt in ('-d', '--distance='):
            distance = int(arg)

    # check wether or not all mandatory parameters have been provided
    if not os.path.isfile(list_of_microRNAs) and os.path.isfile(priority) :
        print_help()
        print "key missing"
        print "Arguments given :", list_of_microRNAs, priority
        sys.exit(1)
        

            
    return



def circular_construction(priorities, set_of_mir_sequence):
    
    #this function initialise the various cluster of binding site and then elongate them iteratively
    
    global current_size
    global list_of_cluster
    global override

    
    for x, priority, mir_seq in zip(xrange(0, len(list_of_cluster)), priorities, set_of_mir_sequence) :
   
       
        for iteration in xrange(0, priority):                       #elongation of a cluster of seed by the number of priority dedicated
            
            sequence=""
            sequence=elongation(sequence, mir_seq)
            current_size+=len(sequence)
            
            if current_size>size :
                
                return

            else :

                if len(override)>=1 and len(override[x])>=10:
                    
                    list_of_cluster[x]+=override[x]
                
                else :
                    
                    list_of_cluster[x]=elongation(list_of_cluster[x], mir_seq)
                
    
    circular_construction(priorities, set_of_mir_sequence)
  

    return
    





if __name__ == '__main__':
    
    reload(sys)
    get_cli(sys.argv[1:])

    filin = open(list_of_microRNAs)
    set_of_mir=set()
    lines=filin.readlines()
    
    for line in lines :
        line=line.replace("\n", "")
        line=line.replace(" ", "")
        set_of_mir.add(line)
    
    
    filin.close()
    
    filin = open(priority)
    priorities=[]
    lines=filin.readlines()
    
    for line in lines:
        line=line.replace("\n", "")
        line=line.replace(" ", "")
        priorities.append(int(line))
        

    if not len(set_of_mir)==len(priorities):
        print "Not the same number of microRNAs and priority"
        print "Size set of microRNAs : ", len(set_of_mir), " Number of priority :", len(priorities)
        sys.exit(1)

    if len(list_of_sequence)>=1:            #in case of override for quality control testing
        filin = open(list_of_sequence)
        lines=filin.readlines()
        for line in lines:
            line=line.replace("\n", "")
            line=line.replace(" ", "")
            override.append(line)
    
        if not len(override)==len(priorities):

            while len(override)!=len(priorities):

                line=""
                override.append(line)

        
    list_of_cluster=[]
    set_of_mir_sequence= []
    
    for mir in set_of_mir:
        #a cluster is here a set of seed for the same micro-RNAs
        #the reason for this construction is that it have been proved that the micro-RNAs have a tendancy to "jump" from a binding site to the following one
        #this, of course, only happen if the binding site are for the same micro-RNAs
        #see A Dynamic Search Process Underlies MicroRNA Targeting, D. Chandradoss et al, 2015, http://dx.doi.org/10.1016/j.cell.2015.06.032
        
        
        cluster="" 
        list_of_cluster.append(cluster)
        
        mir_seq=micro_RNA_sequence_request(mir)
        set_of_mir_sequence.append(mir_seq)    
            
     
    latex_init()
    

    current_size=0
    circular_construction(priorities, set_of_mir_sequence)
    
    
     
        
    for x in xrange(0, len(list_of_cluster)):
        list_of_cluster[x]=poly_U_liker(list_of_cluster[x])
        print list_of_cluster[x] 
    

    latex_output_cluster(list_of_cluster, set_of_mir, set_of_mir_sequence)
    
    for x in xrange(0, len(list_of_cluster)):
        list_of_cluster[x]=from_5_prime_to_3_prime_correction(list_of_cluster[x])

        


    r=len(list_of_cluster)

    iterable=list(xrange(0, r))

    
    
    combinations=itertools.combinations(iterable, r)
    list_of_circular_RNAs=[]

    for combination in combinations:
        sequence=""
        for indices in combination :
            sequence+=list_of_cluster[indices]
        list_of_circular_RNAs.append(sequence)
    
    fichier=open("output_for_vienna.txt", "w")
    for x, circ in zip(xrange(0, len(list_of_circular_RNAs)) ,list_of_circular_RNAs):
        fichier.write(">" + str(x) + "\n")
        fichier.write(circ)

    fichier.close()

    ViennaRNA_RNAfold("output_for_vienna.txt")
    best_circ, best_score, best_seq = RNAfold_grep()
    latex_elong()
    RNAfold_latex_output(best_circ, best_score, best_seq)
    
    mature_reader("mature.fa")      #output human_microRNAs.fasta
    
    
    
    latex_output_miranda_header()
    
    for x, mir in zip(xrange(0, len(list_of_cluster)), set_of_mir ):
        #print "x", x
        sequence=list_of_cluster[x]
        fichier=open("output_for_miranda.txt", "w")
        fichier.write(">cluster_of_seed \n")
        fichier.write(sequence)
        fichier.close()
        
        x="miranda human_microRNAs.fasta output_for_miranda.txt -noenergy -strict -sc 150 > output_miranda"
        os.system(x)     
        
        list_of_dictionnary=miranda_score_reader("output_miranda")
        best_alignement=selecting_best_hit_score(list_of_dictionnary, number_of_best_hits=3)
        latex_output_miranda(best_alignement, mir)
        #print best_alignement

    
    latex_report_pdf()























