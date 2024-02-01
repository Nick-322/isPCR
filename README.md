# isPCR
Package for in silico PCR


This python package performs the following: 

1. Performs isPCR on two provided assembly files with a single provided primer file.
2. Generates one amplicon from each assembly file that is aligned to each other.
3. Checks which orientation the two sequences align best in and only return the best 
     alignment.
4. Prints the alignment followed by the alignment score to terminal.


Sample usage: 

    ./amplicon_align.py -1 data/Pseudomonas_aeruginosa_PAO1.fna -2 data/Pseudomonas_protegens_CHA0.fna -p
    data/rpoD.fna -m 2000 --match 1 --mismatch=-1 --gap=-1
