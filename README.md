# IGACAT: 
## Purpose

Given the 4 text files phaseI_midyearpulse_nt_to_aa_complete.txt, phaseI_midyearpulse_selected_seqs_key_final.txt,
phaseI_midyearpulse_selected_seqs_randomized.fasta, phaseI_midyearpulse_whole_nt.fasta, and a directory that
store subs matrix that we want to use for alignment, it will report all the possible alignment (including 6 frames translation)

## Requirements
* [Python](https://www.python.org/)
* [Biopython](https://biopython.org/wiki/Download)
* [fasta6](https://github.com/wrpearson/fasta36)
## Installation
Users can either use github interface Download button or type the following command in command line:
```bash
git clone https://github.com/nguyenngochuy91/IGACAT.git
```
The users can either download the source codes of the requirements or use package management systems such as [brew](https://brew.sh/),
 or [conda](https://conda.io/miniconda.html). Users can also try [pip](https://pypi.org/project/pip/).

If user has conda, to install [fasta6](https://github.com/wrpearson/fasta36), user can type in the command line:
```bash
conda install -c biobuilds fasta
```


## Usage
The easiest way to run the project is to execute the script [solve.py](https://github.com/nguyenngochuy91/IGACAT/blob/master/solve.py)
One can either run it by typing the following in command line:
```bash
./solve.py
```
or :
```bash
python solve.py
```
The script was quickly written to serve the purpose of getting result,
I will reformat and structure it soon.  
The alignment part and the part to pick out the best result were commented out for time consumming reason.
You are more welcome to run those part. The alignment parts start from line 194-208. The analizing parts start start from 218-221.
What the current script does is to read from file seqToProteinStandard.txt (result of using only the standard translation table)
 and seqToProtein.txt(result using all translation table), and compare them with the key.

seqToProteinStandard.txt file is a dictionary, the key is the fragment sequence (permuted), the value is an array of this format 
(originalProtein,percentIdentity,evalue,substitutionMatrix,TranslationTable). 

The fragments directory store all the translation of our sequence fragment using 43 translation table, in 6 reading frames. 
Each of the file is in fasta format, the id for example could be **Seq5405_166_bp_1_2_1**. Here, **Seq5405_166_bp** indicate the frag seq id,
**1_2_1** indicates this is on normal strand, reading frame 3, and the **second** index in the list. Another one from the same sequence **Seq5405_166_bp_1_2_2**
means almost the same thing except **1_2_2** indicates this is on normal strand, reading frame 3, and the **third** index in the list.
The reason for this is there might be stop codon in the middle of the fragment when translating.

I also generate an **analyze.csv** file, there are 4 tables. 
1. Precions and Recal table: This table has 7 columns
    * Column **protein** stores the protein ID.
    * Column **precision** stores the precision value of the protein ID.
    * Column **recall** stores the recall value of the protein ID.
    * Column **f1-score** stores the f1-score of the protein ID.
    * Column **support** stores the support value (number of permuted fragement mapped to this protein ID).
    * Column **bin** stores the bin threat of the protein ID.
    * Column **name** stores the protein name of the protein ID.
2. Protein table: This table has 2 columns. 
    * Column **protein** stores the protein ID.
    * Column **count** stores the number of time such protein ID
3. Translation table: This table has 2 columns. 
    * Column **translationTable** stores the name of translation table.
    * Column **count** stores the number of time the translation table was used correctly.
4. Substitution table: This table has 2 columns. 
    * Column **subsTable** stores the name of substitution table.
    * Column **count** stores the number of time the substitution table was used correctly.     