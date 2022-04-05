# Transcription-Translation-Reverse_Translation
This repository includes codes for transcription-translation &amp; reverse translation.

## Program Descriptions:
### 1) Transcription & Translation:
The following program prompts the user to input a proper fasta file which contains DNA sequences for several proteins.
The program then reads and scans this file. Then, it uses the stored data to perform several actions like extracting
the number of exons, introns, %AGCT etc... and writes it to a new file called DNAstats.txt alongside displaying it
on the screen for the user. After that, the program prompts the user for a protein ID and checks if that ID is found
within the given sequences. If not it prompts the user again, and if it found a match it then displays all the 
statistics for that protein's DNA and then displays the mRNA sequecne for that specific protein ID and also displays
the sequence of amino acids that this mRNA generates (respecting start & stop codons). Finally, the program takes the
upstream nucleotides (before first exon), reads it in windows of 20 nucleotides and determines the melting temperatures
of the primers, gets the average - range - #of windows - min - and max.
The program terminates right after displaying the execution time.

### 2) Reverse Translation:
This progran takes a peptide sequence as a user input and then displays all the mRNA sequences that can be formed 
from this peptide sequence along with all the %GC of each mRNA from lowest % to highest %. Then, the program finds
the mRNAs that have a %GC closest to 50% and displays them as the most probable mRNAs.
The program terminates after displaying the execution time.
