#!/usr/bin/env python
# coding: utf-8

# # Reverse Translation
# Frederick Abi Chahine

# In[1]:


#!/usr/bin/env python

#Final version of part 2

#Program Description:
#This progran takes a peptide sequence as a user input and then displays all the mRNA sequences that can be formed 
#from this peptide sequence along with all the %GC of each mRNA from lowest % to highest %. Then, the program finds
#the mRNAs that have a %GC closest to 50% and displays them as the most probable mRNAs
#the program terminates after displaying the execution time.

import re                     #regex will be used to ensure that the user input a proper peptide sequence.
from datetime import datetime #will be used in order to display the run time

def generateAllPossibilities(pointer, matrix, mrna_string, possibilities_list):
    #This function is here in order to generate every possible unique mrna sequence from the combination of
    #the codons from the matrix created in function reverseTranscribe. It respects the order of the matrix,
    #in order to generate the appropriate mRNAs for the given peptide sequence.
    
    if pointer==len(matrix):
        possibilities_list.append(mrna_string)
    
    else:
        for codon in matrix[pointer]:
            new_pointer= pointer + 1
            new_mrna_string= mrna_string + codon
            generateAllPossibilities(new_pointer, matrix, new_mrna_string, possibilities_list)
    
    return possibilities_list

def reverseTranscribe(pro_seq):
    #This function takes the protein sequence that the user input as an argument and then generates a matrix in
    #which the mother list has sub-lists, and every sub-list has all possible codons for a specific amino acid.
    #It then computes the # of appropriate mrna sequences that can be generated from this sequence, and then 
    #returns the generated matrix.
    
    codon_map = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
                 "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",          # UCC is lower case 's'. Values may differ if changed to upper case 'S'
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
    
    matrix= []
    
    for i in range(len(pro_seq)):
        matrix.append([key for key,val in codon_map.items() if val == pro_seq[i]])
    
    return matrix #every sub-list has all codons for a specific amino acid

def computePercentGC(all_mrnas):
    #This method take the previously created list that has all the mRNA sequences in it, and creates a new list that
    #stores all the %GC for each mRNA respectively (same indices). Then, with the help of zip() we are able to sort 
    #the %GC list and in turn shadow sort the mRNA list with respect to the indices => every mRNA will be at the 
    #same index of its %GC. We need to sort in order to display the mRNAs with their %GC from lowest to highest.
    #Finally, we return both of the lists to later on display them on the screen for the user.
    
    percent_GC= [] #this list will store all the %GC of every mRNA
    for mrna in all_mrnas: 
        #loops through the mRNA list and calculates the %GC for each mRNA to store it in percent_GC
        gc_per= round(((mrna.count('G')+mrna.count('C'))/len(mrna))*100, 2)
        percent_GC.append(gc_per)
    
    #print(percent_GC)
    #print(all_mrnas)
    
    zipped_lists = zip(percent_GC, all_mrnas)  #first is what it sorts (the %) second goes with it (the mRNA)
    sorted_pairs = sorted(zipped_lists)        #shadow sorting the mRNA while sorting the %GC
    tuples = zip(*sorted_pairs)
    percent_GC, all_mrnas = [ list(tuple) for tuple in  tuples]
    
    return all_mrnas, percent_GC
    
def main():
    
    pro_seq= input("Enter protein sequence: ") #ACDEFGHIKLMNPQRSTVWY   ONLYYYY

    while re.findall(r'[^AC-IK-NP-TVWY]', pro_seq): #To ensure that the user enters any correct combination of amino acids, and avoids any typos/wrong letters/wrong case of letters
        print("\n**You entered an invalid letter/character in your sequence. Only input \"ACDEFGHIKLMNPQRSTVWYs\".**\n")
        pro_seq= input("Enter protein sequence:")
    
    start=datetime.now()                                     #this stores the current time in order to subtract it from the time at the end of execution => we get run time
    matrix = reverseTranscribe(pro_seq)                      #invokes the method in order to create the 2D list
    all_mrnas= generateAllPossibilities(0, matrix, "", [])   #invokes the recursive method to generate all combinations of mRNA
    all_mrnas, percent_GC = computePercentGC(all_mrnas)      #invokes the method that computes and sorts the mrnas with their respective %GC
    print()                                                  #to improve display
    
    for i in range(len(percent_GC)):                         #this loop simply displays the content to the screen
        print(all_mrnas[i], ": %GC =", percent_GC[i])
    
    print("\n---> All possible combinations:", len(all_mrnas), "\n")  #displays the number of mRNA combinations
    print("---> The most probable mRNA sequence(s) is/are: (The ones with %GC closest to 50%)\n")
                
    index_array= []                                    #this array will store all the indices of the mRNAs with %GC closest to 50 (or 50 itself)
    result= abs(percent_GC[0]-50)                      #this will be an initializer for the if condition inside the loop (we get absolute value since we do not want the sign to affect our code)
    for i in range(len(percent_GC)):
        if abs(percent_GC[i]-50) < result:             #the percent - 50 will give us the number it is away from 50, so if the number is smaller, it means that it is closer to 50, if it is 0 then it is 50
            result= abs(percent_GC[i]-50)
            index_array= []                            #we have to refresh the index array if it enters this statement in order to remove all previous indices that were assumed closest to 50
            index_array.append(i)
        elif abs(percent_GC[i]-50) == result:          #else if the distance is equal to the smallest then we simply append it to the array to display later on
            index_array.append(i)
            
    for index in index_array:                          #this loop simply displays onto the screen all the mRNA sequences with %GC closest to 50 (or 50 itself)
        print(all_mrnas[index], ": %GC =", percent_GC[index])
    
    print('\n( Run Time: ', datetime.now()-start, ")") #this displays the execution time of this code

main()


# In[ ]:




