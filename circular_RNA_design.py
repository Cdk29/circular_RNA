#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 13:24:04 2017

@author: erolland
"""

import requests
from StringIO import StringIO
from lxml import etree
from bs4 import BeautifulSoup

def poly_U_liker(sequence):
    #micro-RNA binding site with hight local content of  A and U perform the best 
    #MicroRNAs: Target Recognition and Regulatory Functions
    #DOI 10.1016/j.cell.2009.01.002
    #there is some protein which seems to target polyU so, every 3 U a A is put instead of a U
    
    count=0
    returned_sequence=""
    for i in sequence:
        if i=="X" or i=="W":
            if count==3:
                returned_sequence=returned_sequence+"A"
                count=0
            else:
                returned_sequence=returned_sequence+"U"
                count=count+1
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

def seed_site_distance(micro_RNA_sequence1, sequence):
    # Saetrom et al,(2007) Distance constraints between microRNA target sites dictate efficacy and cooperativity. 
    # DOI 10.1093/nar/gkm133
    # creation of the filling sequence of 13-35 nt between the two microRNA binding site for optimal cooperation
    seed_site_distance=len(micro_RNA_sequence1[8:len(micro_RNA_sequence1)])
    if seed_site_distance < 13 :
        while seed_site_distance != 13:
            sequence=sequence+"X"
            seed_site_distance=seed_site_distance+1
    return sequence



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

def length_control(max_circular_rna_length, sequence, micro_RNA_sequence):
    sequence=seed_site_distance(micro_RNA_sequence, sequence)
    if len(sequence)+len(micro_RNA_sequence)+12 >= max_circular_rna_length:
        return SEQ_GE
    else:
        return 0

def elongation(sequence, micro_RNA_sequence):
    
    sequence = micro_RNA_binding_site(micro_RNA_sequence, sequence)
    sequence = seed_site_distance(micro_RNA_sequence, sequence)
    
    return sequence

def half_length_control(max_circular_rna_length, sequence):
    half=max_circular_rna_length/2
    if len(sequence)>= half :
        return REACH
    else:
        return 0
    
def structure_relax(sequence, max_circular_rna_length):
    while len(sequence) !=  max_circular_rna_length:
        sequence=sequence+"W"
    return sequence


def main_circular_RNA(micro_RNA_1, micro_RNA_2, max_circular_rna_length):
    
    Half=False
    sequence = ""
    micro_RNA_sequence1=micro_RNA_sequence_request(micro_RNA_1)
    micro_RNA_sequence2=micro_RNA_sequence_request(micro_RNA_2)

    while(1):
        sequence=elongation(sequence, micro_RNA_sequence1)
        control=length_control(max_circular_rna_length, sequence, micro_RNA_sequence2)
        if control==SEQ_GE:
            break
        if Half==False:
            control=half_length_control(max_circular_rna_length, sequence)
            if control==REACH:
                sequence=sequence+"W"*12
                Half=True
        sequence=elongation(sequence, micro_RNA_sequence2)
        control=length_control(max_circular_rna_length, sequence, micro_RNA_sequence1)
        if control==SEQ_GE:
            break

    sequence=structure_relax(sequence, max_circular_rna_length)
    sequence=poly_U_liker(sequence)
    sequence=from_5_prime_to_3_prime_correction(sequence)

    return sequence

SEQ_GE = ">="
REACH  = "half of the size of the circular RNA is reach"


micro_RNA_1="hsa-miR-210"
micro_RNA_2="hsa-miR-210"
max_circular_rna_length=900


sequence=main_circular_RNA(micro_RNA_1, micro_RNA_2, max_circular_rna_length)