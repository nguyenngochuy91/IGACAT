#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  4 12:52:18 2018

@author: huyn
"""

from Igacat import Keys, Permutations, Genes
from Bio import Seq
import argparse
from Bio.Data import CodonTable
import os
import json
import time
#this function will return all of the files that are in a directory. os.walk is recursive traversal.
def returnRecursiveDirFiles(root_dir):
    result = []
    for path, dir_name, flist in os.walk(root_dir):
        for f in flist:
            fname = os.path.join(path, f)
            if os.path.isfile(fname):
                result.append(fname)
    return result

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

#     
def getInfo():      
    permutations = Permutations()
    # translate all, this output has the followign format {translation table name : 
    #{fragment name: {strand: [0:{},1:{},2:{}],reverse: [0:{},1:{},2:{}]} }}
#    print ("Translating ...")
    translation       = translateFragments(permutations)
    # genes has record attribute, gene.record has format : {protein id: protein seq}
    genes        = Genes()
    # list of matrix to use
    return translation,genes
# take output from getInfo above
def writeTranslation(translation,genes):
    alignments ="/home/huyn/IGACAT/alignments/"
    for table in translation:
        dir = "/home/huyn/IGACAT/fragments/"+table.replace(" ","_")+"/"
        try:
            os.mkdir(dir)
            os.mkdir(alignments+table+"/")
        except:
            print ("it was already created!!!")
        for fragment in translation[table]:
            outfile = open(dir+fragment+".fasta","w")
            for strand,info in translation[table][fragment].iteritems():
                for orf in range(3):
                    for i in range(len(info[orf])):
                        name = ">{}_{}_{}_{}\n".format(fragment,strand,orf,i)
                        outfile.write(name)
                        sequence = "{}\n".format(info[orf][i])
                        outfile.write(sequence)
            outfile.close()
    outfile = open("genes.fasta","w")
    for g in genes.record:
        outfile.write(">{}\n{}\n".format(g[0],g[1]))
    outfile.close()
## get the bin threat from the file
def getThreatBin():
    proteinThreat = Keys("phaseI_midyearpulse_nt_to_aa_complete.txt")
    binsThreat    = {}
    for item in proteinThreat.record:
        id        = item['SeqID'].split(".")[0]
        threatBin = item['Threat Bin']
        binsThreat[id] = threatBin
    return binsThreat
## from the alignment, get only the one translate from standard table
def getStandard():
    alignments               = returnRecursiveDirFiles("alignments")
    dictionary = {}

    for path in alignments:
        if "Standard" not in path:
            continue
        output = getBest(path)
        if not output:

            continue
        info = path.split("/")[-1]
        seqName = info.split("_")[0]
        
        if seqName in dictionary:
            currentInfo = dictionary[seqName]
            if output[1]>= currentInfo[1] and output[2]<=currentInfo[2]:
                dictionary[seqName] = output
        else:

            dictionary[seqName] = output
    json.dump(dictionary,open("seqToProteinStandard.txt","w"),indent=4)

    return dictionary 

## from all the alignment get the best
def getBest(path):
    info1 = path.split("/")[-1]
    info2 = info1.split("_")
    matrix  = info2[-1]
    translationTable = path.split("/")[-2]
    infile = open(path,"r")
    currentMin = 10000
    currentMax = 0
    output = None
    for line in infile.readlines():
        line = line.split('\t')
        protein = line[1].split('|')[-1].split(".")[0]
        percentId = float(line[2])
        eVal      = float(line[-2])
        if eVal<=currentMin and percentId>=currentMax: 
            output = (protein,percentId,eVal,matrix,translationTable)
            currentMin = eval
            currentMax = percentId
    infile.close()
    return output

def getDictionary():
    alignments               = returnRecursiveDirFiles("alignments")
    dictionary = {}

    for path in alignments:
        print (path)
        output = getBest(path)
        if not output:

            continue
        info = path.split("/")[-1]
        seqName = info.split("_")[0]
        
        if seqName in dictionary:
            currentInfo = dictionary[seqName]
            if output[1]>= currentInfo[1] and output[2]<=currentInfo[2]:
                dictionary[seqName] = output
        else:

            dictionary[seqName] = output
    json.dump(dictionary,open("seqToProtein.txt","w"),indent=4)

    return dictionary 


def compare(seqToProtein,reduceKey):
    count = 0
    missing = []
    for k in seqToProtein:
#        print (reduceKey[k],seqToProtein[k])
        if reduceKey[k] in seqToProtein[k][0]:
            count+=1
        else:
            missing.append([k,reduceKey[k],seqToProtein[k][2]])
    return float(count)/len(seqToProtein),missing
def main():
    start = time.time()
    translation,genes = getInfo()
    writeTranslation(translation,genes)
    substitutionMatrixFiles = returnRecursiveDirFiles("substitutions")
    fragments               = returnRecursiveDirFiles("fragments")
    genePath                = "genes.fasta"
    # generate the alignment 
#    for pathToFrag in fragments:
#        infoToFrag = pathToFrag.split("/")
#        for pathSubMatrix in substitutionMatrixFiles:
#            matrixName       = pathSubMatrix.split("/")[-1]
#            fragmentName     = infoToFrag[-1].replace(".fasta","")
#            translatioNTable = infoToFrag[1]
#            dir  = "/home/huyn/IGACAT/alignments/{}/".format(translatioNTable.replace(" ","_")+"/")
#            try:
#                os.mkdir(dir)
#            except:
#                pass
#            name = "{}{}_{}".format(dir,fragmentName,matrixName)
#            cmd  = "ssearch36 -T 8 -m 8 -s {} {} {} > {}".format(pathSubMatrix,pathToFrag,genePath,name)
#            print ("Working on {} ...".format(cmd))
#            os.system(cmd)
#     for each protein in they key, get the threat bin of it
#    binsThreat = getThreatBin()
#    # from our alignment of different translation table, and different subs matrix, store the best in seqToProteinInfo.txt
#    dictionary = getDictionary()
#    # from the above, map seq to Threat
##    seqToThreat = mapSeqThreat(binsThreat,dictionary)
#    # standard table
#    standard    = getStandard()
    seqToProtein = json.load(open("seqToProtein.txt","r"))
    seqToProteinStandard = json.load(open("seqToProteinStandard.txt","r"))
    # getting the key result, to make comparison
    k = Keys("phaseI_midyearpulse_selected_seqs_key_final.txt")
    # what we really interested only mapping fragment seq to original gene
    reduceKey = {}
    for item in k.record:
        reduceKey[item["SeqID"]] = item['Gene']
    # compare between only using the standard translation and the key
    accuracyStandard,missingStandard = compare(seqToProteinStandard,reduceKey)
    print ("accuracyStandard:",accuracyStandard)
    # compare between only using all translation table and the key
    accuracyAll,missing = compare(seqToProtein,reduceKey)
    print ("accuracyAll:",accuracyAll)
    stop = time.time()
    print ("Wow!, it took us {} secons".format(stop-start))
main()