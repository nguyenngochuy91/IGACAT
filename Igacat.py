#!/usr/bin/env python
''' Author  : Huy Nguyen
    Program : IGACAT challenge
              
    Start   : 05/08/2016
    End     : /08/2016
'''
from Bio import SeqIO
from Bio import Entrez
import time
class Permutations(object):
    def __init__(self,fasta= "phaseI_midyearpulse_selected_seqs_randomized.fasta"):
        try:
            self.record = list(SeqIO.parse(fasta,"fasta"))
        except:
            self.record = []
            print ("Either file is mal-formated, or different version of Python\n")
            
class Genes(object):
    def __init__(self,fasta= "phaseI_midyearpulse_whole_nt.fasta",flag= False):
        try:
            record = list(SeqIO.parse(fasta,"fasta"))
            self.record = []
            for r in record:
                seq = r.seq
                if not flag:
                    seq = seq.translate().strip("*")
                    seq = str(seq).replace("*","")
                id  = r.id
                self.record.append([id,seq])
        except:
            self.record = []
            print ("Either file is mal-formated, or different version of Python\n")

class Keys(object):
    def __init__(self,fasta= "phaseI_midyearpulse_whole_nt.fasta" ):
        self.record,self.title = self.parse(fasta)
    def parse(self,fasta):
        record = []
        handle = open(fasta,"r")
        title = handle.readline().strip().split('\t')
        for l in handle.readlines():
            l    = l.strip()
            d    = {}
            info = l.split('\t')
            for i in range(len(title)):
                d[title[i]] = info[i]
            record.append(d)
        return record,title

    
    
#k = Keys(fasta= "phaseI_midyearpulse_selected_seqs_key_final.txt")
#gene = set([item["Gene"] for item in k.record])
#outfile = open("phaseI_midyearpulse_whole_aa.fasta","w")
#for g in gene:
#    # create a handle to fetch from ncbi
#    print (g)
#    handle = Entrez.efetch(db="nucleotide", id=g, rettype="fasta", retmode="text")
#    lines  = handle.read()
#    outfile.write(lines.strip())
#    outfile.write("\n")
#    time.sleep(.4)
#outfile.close()