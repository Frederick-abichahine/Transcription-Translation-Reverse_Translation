#!/usr/bin/env python
# coding: utf-8

# # Transcription - Translation
# Frederick Abi Chahine

# In[28]:


#!/usr/bin/env python

#Final version of part 1

#Program Description:
#The following program prompts the user to input a proper fasta file which contains DNA sequences for several proteins.
#The program then reads and scans this file. Then, it uses the stored data to perform several actions like extracting
#the number of exons, introns, %AGCT etc...and writes it to a new file called DNAstats.txt alongside displaying it
#on the screen for the user. After that, the program prompts the user for a protein ID and checks if that ID is found
#within the given sequences. If not it prompts the user again, and if it found a match it then displays all the 
#statistics for that protein's DNA and then displays the mRNA sequecne for that specific protein ID and also displays
#the sequence of amino acids that this mRNA generates (respecting start & stop codons). Finally, the program takes the
#upstream nucleotides (before first exon), reads it in windows of 20 nucleotides and determines the melting temperatures
#of the primers, gets the average - range - #of windows - min - and max.
#The program terminates right after displaying the execution time.

import re                               #we import regular expression in order to use it to find all lines that start with ">" since they indicate a start to a new sequence
import os                               #used to ensure proper file input from user
from datetime import datetime           #will be used in order to display the run time

count_seq=0                             #a counter to count and display the number of sequences in file
length= 0                               #a variable used to increment the lengths of exons/introns
id_array= []                            #this array will store all the sequence IDs from the file
num_exons= []                           #this array stores the number of exons in each sequence respectively
num_introns= []                         #this array stores the number of introns in each sequence respectively
avg_exon_len= []                        #this array stores the average length of exons in each sequence respectively
avg_intron_len= []                      #this array stores the average length of introns in each sequence respectively
exon_string= ""                         #this is a temporary string to add each exon to it, will be refreshed each loop
exon_array= []                          #stores all exons of each sequence
intron_string= ""                       #this is a temporary string to add each intron to it, will be refreshed each loop
upstream_downstream_intron_string= ""   #similar to the temp exon string, but for all lower case letters (including upstream, downstream, and introns)
upstream_downstream_intron_array= []    #first element of array will be the upstream lower case letters of the sequence; last element of the array will be the downstream lower case letters of the sequence, and all elements in between will be the introns... 
all_dict= {}                            #this dictionary stores all sequence IDs as keys, and their values are the exon, intron, upstream, and downstream nucleotides respectively

files_in_dir=os.listdir()
file_name= input("Enter Name of File: ")

while (file_name not in files_in_dir) or (os.path.isdir(file_name)) or (not file_name.endswith('.fasta')):
    #This loop ensures that the user inputs the correct file format.
    #a corrcet file should be a 1) .fasta file 2) present in working directory 3) a readable file/not a directory
    if file_name not in files_in_dir:
        print("\n**File missing / not in directory. Try again**\n")
    elif os.path.isdir(file_name):
        print("\n**File is a directory. Try again**\n")
    elif not file_name.endswith('.fasta'):
        print("\n**File needs to have an extension \".fasta\". Try again**\n")
    file_name= input("Enter Name of File: ")

start=datetime.now()          #this stores the current time in order to subtract it from the time at the end of execution => we get run time
with open(file_name) as file: #"FinalProject.fasta"
    for line in file:
        if re.findall(r'>', line):                                                 #checks if the lines start with ">"
            if (len(exon_array)>0 and len(upstream_downstream_intron_array)>0):    #we only want to displays and refresh the arrays at the END of each sequence. However, the last sequence wont enter this condition => we account for it after the for loop!
                upstream_downstream_intron_array.append(upstream_downstream_intron_string)
                upstream_downstream_intron_string= ""
                num_exons.append(str(len(exon_array)))                             #will give us the number of exons in this sequence
                num_introns.append(str((len(upstream_downstream_intron_array)-2))) # -2 since we do not want to include the upstream and downstream sequences
                
                for exon in range(len(exon_array)):                                #this loop will give us the total sum of the length of exons
                    exon_string+= exon_array[exon]
                    length+= len(exon_array[exon])
                    
                avg_exon_len.append(str(round(length/len(exon_array))))            #calculates the average length of the exons
                length= 0                                                          #refreshes length to be used for introns
                
                for intron in range(1, len(upstream_downstream_intron_array)-1):   #this loop will give us the total sum of the length of introns, we start from 1 to skip index 0 which is the upstream sequence and we range to len-1 to stop before the downstream sequence
                    intron_string+= upstream_downstream_intron_array[intron]
                    length+= len(upstream_downstream_intron_array[intron])
                    
                avg_intron_len.append(str(round(length/(len(upstream_downstream_intron_array)-2))))
                all_dict[id_array[count_seq-1]]= exon_string, intron_string, upstream_downstream_intron_array[0], upstream_downstream_intron_array[len(upstream_downstream_intron_array)-1]
                #The below lines are to refresh the strings, arrays, and counters
                exon_string=""
                intron_string=""
                upstream_downstream_intron_array= []
                exon_array= []
                length= 0
                
            id_array.append(line[1:].strip()) #adds the ID of the sequence to the respective array without the leading '>'
            count_seq+=1                      #if they do start with ">" then we increment the counter
            
        else:                                                      #if the line is not the sequence ID it enters here
            for nt in line:                                        #checks every letter of that given line
                if (nt=='a' or nt=='g' or nt=='c' or nt=='t'):     #if the letter is lower case, it should be added to its respective string and appended later on to the right array
                    upstream_downstream_intron_string+= nt
                    
                    if len(exon_string) > 0:                       #if it is not the first loop and the exon string is not empty, we know that we just finished reading an exon, hence we must append to the exon array and refresh the exon string for the next exon (if any) in this sequence
                        exon_array.append(exon_string)
                        exon_string= ""                            #refresh the string
                        
                elif (nt=='A' or nt=='G' or nt=='C' or nt=='T'):   #if the letter is upper case, it should be added to its respective string and appended later on to the right array
                    exon_string+= nt
                    
                    if len(upstream_downstream_intron_string) > 0: #if the string for lower case letters is not empty then we must append that string to its respective array and refresh it as we know that we have just finished reading it from the file
                        upstream_downstream_intron_array.append(upstream_downstream_intron_string)
                        upstream_downstream_intron_string= ""      #refresh the string
    
    #since the loop does not account for the very last sequence, we shall append them here
    
    upstream_downstream_intron_array.append(upstream_downstream_intron_string)
    upstream_downstream_intron_string= ""
    num_exons.append(str(len(exon_array)))
    num_introns.append(str((len(upstream_downstream_intron_array)-2)))
    
    for exon in range(len(exon_array)):
        exon_string+= exon_array[exon]
        length+= len(exon_array[exon])
        
    avg_exon_len.append(str(round(length/len(exon_array))))
    length= 0
    
    for intron in range(1, len(upstream_downstream_intron_array)-1):
        intron_string+= upstream_downstream_intron_array[intron]
        length+= len(upstream_downstream_intron_array[intron])
        
    avg_intron_len.append(str(round(length/(len(upstream_downstream_intron_array)-2))))
    all_dict[id_array[count_seq-1]]= exon_string, intron_string, upstream_downstream_intron_array[0], upstream_downstream_intron_array[len(upstream_downstream_intron_array)-1]
    
    print("\n---> The file contains", count_seq, "sequence(s).\n") #just for display currently
    
    index= 0 #This index is crucial for the proper output and appending in the following loop...
    writeFile= open("DNAstats.txt", 'w')
    header= (f'{"SequenceID":<15s}\t{"#Exons":<15s}\t{"#Introns":<15s}\t{"AvgExonLength":<15s}\t{"AvgIntronLength":<15s}\t{"%A In Exon":<15s}\t{"%C In Exon":<15s}\t{"%G In Exon":<15s}\t{"%T In Exon":<15s}\t{"%A In Intron":<15s}\t{"%C In Intron":<15s}\t{"%G In Intron":<15s}\t{"%T In Intron":<15s}\n')
    #header= (f'{"SequenceID":<15s}{"#Exons":<15s}{"#Introns":<15s}{"AvgExonLength":<15s}{"AvgIntronLength":<15s}{"%A In Exon":<15s}{"%C In Exon":<15s}{"%G In Exon":<15s}{"%T In Exon":<15s}{"%A In Intron":<15s}{"%C In Intron":<15s}{"%G In Intron":<15s}{"%T In Intron":<15s}\n')
    
    writeFile.write(header)
    writeFile.close()
    writeFile= open("DNAstats.txt", 'a')
    
    for key in all_dict:
        #This loop calculates the percent of every unique nucleotide inside the exon and intron of every sequence
        #It then appends the previously derived data with the following calculations into a new text file
        exon= all_dict[key][0]     #The values of the dictonary stores exons at index 0
        intron= all_dict[key][1]   #The values of the dictonary stores introns at index 1
        percent_A_exon= str(round(((exon.count('A'))/len(exon))*100, 3))
        percent_C_exon= str(round(((exon.count('C'))/len(exon))*100, 3))
        percent_G_exon= str(round(((exon.count('G'))/len(exon))*100, 3))
        percent_T_exon= str(round(((exon.count('T'))/len(exon))*100, 3))
        percent_a_intron= str(round(((intron.count('a'))/len(intron))*100, 3))
        percent_c_intron= str(round(((intron.count('c'))/len(intron))*100, 3))
        percent_g_intron= str(round(((intron.count('g'))/len(intron))*100, 3))
        percent_t_intron= str(round(((intron.count('t'))/len(intron))*100, 3))
        
        append= (f'{key:<15s}\t{num_exons[index]:<15s}\t{num_introns[index]:<15s}\t{avg_exon_len[index]:<15s}\t{avg_intron_len[index]:<15s}\t{percent_A_exon:<15s}\t{percent_C_exon:<15s}\t{percent_G_exon:<15s}\t{percent_T_exon:<15s}\t{percent_a_intron:<15s}\t{percent_c_intron:<15s}\t{percent_g_intron:<15s}\t{percent_t_intron:<15s}\n')
        #append= (f'{key:<15s}{num_exons[index]:<15s}{num_introns[index]:<15s}{avg_exon_len[index]:<15s}{avg_intron_len[index]:<15s}{percent_A_exon:<15s}{percent_C_exon:<15s}{percent_G_exon:<15s}{percent_T_exon:<15s}{percent_a_intron:<15s}{percent_c_intron:<15s}{percent_g_intron:<15s}{percent_t_intron:<15s}\n')
        writeFile.write(append)
        index+=1
    
    writeFile.close()
    name= "DNAstats.txt"
    readFile= open(name, 'r')
    fileContent= readFile.read()
    print(fileContent,"\n")        #This is to display the content/statistics onto the screen for the user
    readFile.close()
    print("---> File called \"" + name + "\" successfully created with all the statistics!\n")
    
    codon_map = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
                 "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
                 "UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP",
                 "UGU":"C", "UGC":"C", "UGA":"STOP", "UGG":"W",
                 "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
                 "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
                 "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
                 "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
                 "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
                 "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
                 "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
                 "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
                 "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
                 "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
                 "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
                 "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G"}
    
    id_input= input("Enter Sequence ID: ").upper()
    marker= True                                            #This is a flag that will be used to indicate if this ID was found or not...
    
    while (marker):                                         #As long as the ID is not found, marker will remain True and it will keep looping
        with open(name) as file:
            head= file.readline()
            
            for line in file:                    
                line_array= line.split("\t")                #stores that line into an array split by the tab, to check if the ID at index 0 matches what the user input
                
                for i in range(len(line_array)):
                    line_array[i]= line_array[i].strip(' ') #we need to strip any spaces to have an accurate check
                                                            #line_array[0] will hold the ID of the sequence
                if id_input == line_array[0].upper():       #If this evaluates to true then the user input was correct
                    marker= False                           #changes the marker to false in order to indicate that the ID entered by the user was correct and to stop looping in the mother while loop 
                    print("\n" + head.strip())              #displays the header for the user
                    print(line)                             #displays the statistics of that sequence for the user
                    exon=  all_dict[line_array[0]][0]       #we know from before that the dictionary has all the exons of that specific sequecne combined at value index 0
                    mRNA = exon.replace("T", "U")           #here we convert the exons into mRNA by replacing all "T" with "U"
                    print("---> The mRNA for the sequence is:\n\n" + mRNA + "\n")
                    
                    index=0                                 #This index will help us in 2 things: 1) to find the start codon & 2) to keep track of codons which are 3 nucleotides
                    peptide_seq= ""                         #will store the entire translated peptide sequence for this mRNA
                    start_codon= True                       #a marker to make sure a start codon was found, we change to false if no start codon was found
                    for i in range(len(mRNA)-2):            #This loop is simply here to find the first start codon and to increment index until it reaches to that codon
                        check_start_codon= mRNA[i]+mRNA[i+1]+mRNA[i+2]
                        if check_start_codon == "AUG":      #i+1 and i+2 in order to get 3 nucleotides to form a codon
                            break                           #we break out as soon as we find the start codon with the appropriate index
                        else:                              
                            index+=1                        #else if the codon is not a start codon then we increment index and continue the loop
                            if (i==len(mRNA)-3):
                                print("---> No start codon found...")
                                start_codon= False          #this indicates that no start codon was found => we shall no perforn any of the following loops
                            
                    while (start_codon) and ((index+3 < len(mRNA)) and (codon_map[mRNA[index:index+3]] != "STOP")): 
                        #this while loop keeps looping and adding to the peptide_seq until a stop codon is reached or we reach the end of the mRNA
                        peptide_seq+= codon_map[mRNA[index:index+3]] # +3 because we want the proper codon
                        index+=3                                     #we increment index by 3 since codons are a group of 3 nucleotides and we must respect that to get the proper peptide sequence
                    
                    if (start_codon):
                        print("---> The peptide sequence is:\n\n" + peptide_seq)
                    break
                    
            if(marker):                                              #if ID that user input did not match any ID present, then marker will remain true and enter this if statement which prompts the user for another ID 
                print("\n**Sequence ID Not Found. Try again.**")     #it will keep looping and prompting until the user inputs a proper ID
                id_input= input("Enter Sequence ID: ").upper()
                    
    if (len(all_dict[line_array[0]][2]) <= 20):  #This condition is true when the length of the upstream nucleotide sequecne is 2 or less since that means we only have 1 window
        print("\n---> Since upstream nucleotides are <= 20, this means there is only 1 window.")
        count= 1                                 #count set to 1 here since we use count as the measure for windows
        up= all_dict[line_array[0]][2]           #the upstream sequecne is stores as the 3rd value of the dictionary key
        tm= (up.count('a') + up.count('t'))*2 + (up.count('c') + up.count('g'))*4
        tm_min= tm
        tm_max= tm
        tm_range= 0                              #0 since only 1 window (max - min)
        tm_average= tm
    else:                                        #else if the len of the upstream nucleotide sequence is >20 we enter here since we have more than 1 window
        tm_array= []                             #this array will store the tm for every window in order to get the average, range, min, and max
        up= all_dict[line_array[0]][2]           #the upstream sequecne is stores as the 3rd value of the dictionary key
        count=0                                  #this will serve as a counter to count the number of windows
        for i in range(len(up) - 19):            #loops from 0 to the length of the upstream seq - 19, we do -19 since we are taking frames of 20 and we dont want to get index out of bound error
            window= up[i:i+20]                   #we take the window to be the substring from i to i+20 since i+20 is not included the window will be size 20
            tm= (window.count('a') + window.count('t'))*2 + (window.count('c') + window.count('g'))*4
            tm_array.append(tm)
            count+=1                             #increment counter for windows
        tm_min= min(tm_array)                    
        tm_max= max(tm_array)
        tm_range= tm_max-tm_min
        tm_average= round(sum(tm_array)/len(tm_array), 3)
            
    print("\n---> Melting temperature of the primers for this sequence:\n")
    print("- Window(s):", count)
    print("- Range:", tm_range,"°C")
    print("- Min:", tm_min,"°C")
    print("- Max:", tm_max,"°C")
    print("- Average:", tm_average,"°C")
    print('\n\n( Run Time: ', datetime.now()-start, ") \n~Time will vary depending on user input speed.") #this displays the execution time of this code

