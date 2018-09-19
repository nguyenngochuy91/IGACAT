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
 or [conda](https://conda.io/miniconda.html). Users can also try [pip](https://pypi.org/project/pip/)


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
First of, I want to sorry for the unorganizing script. It was quickly written to serve the purpose of getting result.
I will reformat and structure it soon. 
The alignment part and the part to pick out the best result were commented out for time consumming reason.
You are more welcome to run those part. The alignment parts start from line 194-208. The analizing parts start start from 218-221.
What the current script does is to read from file seqToProteinStandard.txt (result of using only the standard translation table)
 and seqToProtein.txt(result using all translation table), and compare them with the key.
 
seqToProteinStandard.txt file is a dictionary, the key is the fragment sequence (permuted), the value is an array of this format 
(originalProtein,percentIdentity,evalue,substitutionMatrix,TranslationTable). 


