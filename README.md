# Circular RNA design

This script is part of a workflow developped at the Systems Biology and Bioinformatics department, Universität Rostock, as an extension of the TriplexRNA database, https://triplexrna.org/, a database of cooperating microRNAs and their mutual target.

The script design a circular microRNAs sponge based on microRNAs names passed as input, and output a pdf report detailling the steps of the design and a justification to it.

# Required installations


- miRanda is used by the script for quality control. The script come with the file mature.fa, a fasta format sequences of all mature miRNA sequences, download from miRbase (http://www.mirbase.org/ftp.shtml), during may 2017. The miRanda alignment program can be download here : http://www.microrna.org/microrna/getDownloads.do

- ViennaRNA is used for minimum free energy calculation, and can be installed from here : https://www.tbi.univie.ac.at/RNA/. While this script is in Python, there is no need for the python wrapper.



# How to use it

The script design a circular RNA to sponge one or several microRNAs which names are passed as input. It takes as input :

- a file  containing the names of the microRNAs that the user want the sponge to decoy (mandatory). The names have to be the same as in the triplexRNA database. This file is a txt file with just the name of the microRNAs, one per line.

- a file containining a priorities for each microRNAs (mandatory). Basically, the script build the sponge in a recursive way, and add as many binding site as the priority number allocated to each microRNA on each recursion. For example, for a sponge with twice the number of binding site for one microRNAs compare to the other, the file will contain the integer 2 and 1. This file is a txt file with just one integer per line, and have to be in the same order than the microRNAs they refer to.

- a distance of cooperativity (optionnal), passed as argument in command line. This distance, as define in Saetrom et al (2007), between microRNAs target sites dictate the cooperativity in the microRNAs binding to the sponge. Default cooperativity distance : 17.

- a size in nucleotide (optionnal), passed as argument in command line. This number set the maximum size in nucleotides that the sponge should no trespassed. Default size : 300.

- a file of list of sequence for override (optionnal). The override procedure allow the user to specify the sequence of a binding site he desire for one microRNAs. This option is meant to be used in case of bad quality control, for example when a another microRNAs target with a better affinity the binding site. Cf pdf report for more details. An example of a file this procedure is provide with the script.

# Example

The circular design script is provided with three example files :  Example_Mir.txt, Example_Priority.txt and Example_Override.

To construct a circular sponge for two microRNAs, with a priority of 2 and 1 respectively, with no additionnal options enable, one can simply run :


	python circular_RNA_design.py -l Example_Mir.txt -p Example_Priority.txt  

or

	python circular_RNA_design.py -list_of_microRNAs=Example_Mir.txt -priority=Example_Priority.txt


To show the help and see the other options :


	python circular_RNA_design.py -h

When specifying a distance of cooperativity (here 30 nucleotide from a seed to another) :


	python circular_RNA_design.py -l Example_Mir.txt -p Example_Priority.txt -d 30


When specifying a maximum size of 500 nucleotide for the circular sponge :


	python circular_RNA_design.py -l Example_Mir.txt -p Example_Priority.txt -s 500


And finaly, an example of the override procedure, to specify a specific binding site to be repeat for the second microRNAs only :


	python circular_RNA_design.py -l Example_Mir.txt -p Example_Priority.txt -q Example_Override.txt



# About

This circular RNA design script has been developped at the Systems Biology and Bioinformatics department - Universität Rostock, by E. Rolland and A. Bagnacani, under the supervision of the Dr S. Gupta and the Pr O. Wolkenhauer.
