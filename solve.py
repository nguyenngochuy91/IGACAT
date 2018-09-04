#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  4 12:52:18 2018

@author: huyn
"""

from Igacat import Keys, Permutations, Genes
from Bio import Seq
from Bio.SubsMat import MatrixInfo
import argparse

def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--permutations","-i", help="Permuted sequence file",default = "phaseI_midyearpulse_selected_seqs_randomized.fasta")
    parser.add_argument("--genes","-g", help="Whole sequence", default = "phaseI_midyearpulse_whole_nt.fasta")
    parser.add_argument("--output","-o", help="Output file", default = "output.txt")
    args = parser.parse_args()
    return args

def translate(record,table):
    output = {1:{},-1:{}}
    for strand, nuc in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
        for frame in range(3):
            length = 3 * ((len(record)-frame) // 3)
            output[1][frame] = []
            for pro in record[frame:frame+length].translate(table).split("*"):
                output[1][frame].append(pro)
    return output
if __name__ == "main":
    args         = get_arguments()
    permutations = args.permutations
    genes        = args.genes
    output       = args.output
    permutations = Permutations(permutations)
    genes        = Genes(genes)
    # list of matrix to use
    blosum60     = MatrixInfo.blosum60
    blosum100    = MatrixInfo.blosum100
    pam30        = MatrixInfo.pam30
    pam60        = MatrixInfo.pam60
    pam90        = MatrixInfo.pam90