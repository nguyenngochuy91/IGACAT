#!/usr/bin/env python
''' Author  : Huy Nguyen
    Program : IGACAT challenge
              
    Start   : 05/08/2016
    End     : /08/2016
'''
from Bio import SeqIO

class Permutations(object):
    def __init__(self,fasta= "phaseI_midyearpulse_selected_seqs_randomized.fasta"):
        try:
            self.record = list(SeqIO.parse(fasta,"fasta"))
        except:
            self.record = []
            print ("Either file is mal-formated, or different version of Python\n")
            
class Genes(object):
    def __init__(self,fasta= "phaseI_midyearpulse_whole_nt.fasta" ):
        try:
            self.record = list(SeqIO.parse(fasta,"fasta"))
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
    
