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
from Bio.Data import CodonTable
import json
from Bio import pairwise2
from multiprocessing import Pool

import time
def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--permutations","-i", help="Permuted sequence file",default = "phaseI_midyearpulse_selected_seqs_randomized.fasta")
    parser.add_argument("--genes","-g", help="Whole sequence", default = "phaseI_midyearpulse_whole_nt.fasta")
    parser.add_argument("--output","-o", help="Output file", default = "output.txt")
    args = parser.parse_args()
    return args
# get translation of 6 different frames
def translate(record,table):
    output = {1:{},-1:{}}
    for strand, nuc in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
        for frame in range(3):
            length = 3 * ((len(record)-frame) // 3)
            output[strand][frame] = []
            for pro in nuc[frame:frame+length].translate(table = table).split("*"):
                sequence = str(pro)
                if len(sequence)>=10:
                    output[strand][frame].append(str(pro))
    return output
###############################################################################
# translate all of out fragment
def translateFragments(permutations):
    # get all translate tables (43 of them)
    translationTables = CodonTable.unambiguous_dna_by_name.keys()
    output       = {}
    for table in translationTables:
        output[table] = {}
        for r in permutations.record:
            name = r.id
            dictionary  = translate(r,table)
            output[table][name]=dictionary
    return output


def multiRunWrapper(args):
    return singleAlign(*args)
def singleAlign(seq1,seq2,matrix,geneName):
    alignment=  pairwise2.align.localds(seq1,seq2,matrix,-10,-1)[0]
    return (geneName,alignment[2])
# given a fragment sequence, for each reading frame, align with all the protein sequence with a given substition matrix
 # fragment {-1: {0: ['PYQIPSLSYICYLPSALTFRYGQPPTVGDESCF'],
#  1: ['GSNIVLQSSIGGVNCTIYGVNLLNPTSSPA', 'RSVMVNPPPWETSDAS'],
#  2: ['ARTSFFNLASEVLIVPFTA', 'TCLTLPDPQPKLYLLSAIGINVPLWSTPHRGSREMLP']},
# 1: {0: ['WKHLSSPTVGGWP'],
#  1: ['GSISRLPRWGVDHNGTLMPMADNSYNLGWGSGSVKQVYAVNGTINTSDASLKNDVRA'],
#  2: ['EASLVSHGGGLTITER', 'AGDLVGLSSFTP', 'MVQLTPPMLDWSTMFEP']}}
def alignAll(matrix, fragmentSeq,genes):
    dictionary = {}
    # this dictionary store info of which strand the fragment we worry about, then which
    # protein and reading frame that gives the highest score
    for strand in fragmentSeq:
        print ("strand {}".format(strand))
        dictionary[strand]= {}
        for i in range(3):
            potential = fragmentSeq[strand][i]
            for seq1 in potential:
                pool = Pool()
                myList = [[seq1,seq2,matrix,geneName] for geneName,seq2 in genes.record]
                result = pool.map(multiRunWrapper,myList)
                best   = max(result, key = lambda x:x[1])
                print (i,best[0],best[1])  
                dictionary[strand]= (i,best[0],best[1])      
    return dictionary

    
def main():
    args         = get_arguments()
    permutations = args.permutations
    genes        = args.genes
    output       = args.output
    permutations = Permutations(permutations)
    # genes has record attribute, gene.record has format : {protein id: protein seq}
    genes        = Genes(genes)
    # list of matrix to use
    substitutionMatrices = [MatrixInfo.blosum60,
                            MatrixInfo.blosum62,
                            MatrixInfo.blosum65,
                            MatrixInfo.blosum70,
                            MatrixInfo.blosum75,
                            MatrixInfo.blosum80,
                            MatrixInfo.blosum85,
                            MatrixInfo.blosum90,
                            MatrixInfo.blosum95,
                            MatrixInfo.feng,
                            MatrixInfo.fitch,
                            MatrixInfo.genetic,
                            MatrixInfo.gonnet,
                            MatrixInfo.grant,
                            MatrixInfo.ident,
                            MatrixInfo.johnson,
                            MatrixInfo.levin,
                            MatrixInfo.mclach,
                            MatrixInfo.miyata,
                            MatrixInfo.nwsgappep,
                            MatrixInfo.pam120,
                            MatrixInfo.pam180,
                            MatrixInfo.pam250,
                            MatrixInfo.pam30,
                            MatrixInfo.pam300,
                            MatrixInfo.pam60,
                            MatrixInfo.pam90,
                            MatrixInfo.rao,
                            MatrixInfo.risler,
                            MatrixInfo.structure]
    names = MatrixInfo.available_matrices[10:]
    # translate all, this output has the followign format {translation table name : 
    #{fragment name: {strand: [0:{},1:{},2:{}],reverse: [0:{},1:{},2:{}]} }}
    print ("Translating ...")
    translation       = translateFragments(permutations)
#    json.dump(output,open("output.txt","w"))
    # aligner function
    dictionary = {}
    # for each substitution matrix
    print ("Aligning ... ")
    start = time.time()
    for i in range(len(substitutionMatrices)):
        substitutionTable  = names[i]
        matrix             = substitutionMatrices[i]
        dictionary[substitutionTable]= {}
        print ("Aligning using substitutionTable ... {}".format(substitutionTable))
        for translationTable in translation:
            print ("with translation table {} ... ".format(translationTable))
            fragmentSeqs = translation[translationTable]
            dictionary[substitutionTable][translationTable] = {}
            for fragmentName,fragmentInfo in fragmentSeqs.iteritems():
                print ("with sequence name :{}".format(fragmentName))
                dictionary[substitutionTable][translationTable][fragmentName]= alignAll(matrix, fragmentInfo,genes)
            stop = time.time()
            print (stop-start)
            return dictionary
dictionary = main()