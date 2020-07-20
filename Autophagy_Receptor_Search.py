'''
This script checks the amino acid sequences of all the ATG8 hits from our
Autophagy BioID dataset and determines if they possess an LIR motif and
ubiquitin-binding domain.
'''

#import modules
import requests
import sys
import re
from collections import defaultdict

# The permutate function creates a list of all possible combinations of a motif


def permutate(bases, results, depth=0, current_result=""):
    if depth == len(bases):
        results.append(current_result)
    else:
        for base in bases[depth]:
            permutate(bases, results, depth + 1, current_result + base)


# Motif definitions
'''
------------------------------------------------------------------------
#LIR motif search
#LIR consensus is [W/F/Y]-X-X-[L/I/V]
------------------------------------------------------------------------
'''
LIRmotif_aa = [
    "WFY",  # options for first position
    # options for second position; Note: this is all possible amino acids
    "ACDEFGHIKLMNPQRSTVWY",
    "ACDEFGHIKLMNPQRSTVWY",    # ...
    "LIV"]
LIR_motifs = []
permutate(LIRmotif_aa, LIR_motifs)

# open file and generate a set of unique proteins
dict1 = defaultdict(list)

with open("ATG8 preys- ID and names.txt") as f1:
    for line in f1:
        if line.strip():
            a, b = line.strip().split()
            dict1[a].append(b)

f1_set = set(dict1.keys())

# Creating a new text file for outputting the results
my_file = open("Analysis of ATG8 BioID Preys.txt", "w")
my_file.write(
    "Identification of putative LIR motifs and ubiquitin-binding in the ATG8 preys\n\n")

# retrieve each protein sequence and analyze GO annotation and aa sequence
for f in f1_set:
    # check to see if protein binds ubiquitin and if yes, continue to LIR
    # search
    requestURL = "https://www.uniprot.org/uniprot/?query=" + \
        str(f) + "+AND+organism:9606&sort=score&limit=1&columns=go(molecular function)&format=tab"
    r = requests.get(requestURL, headers={"Accept": "text/plain"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    GO = r.text.rstrip("\n")
    if "ubiquitin binding" in GO:
        my_file.write(
            str(set(dict1[f])) + " is a ubiquitin-binding protein based on GO annotation\n")
        # Search for LIR motif
        requestURL = "https://www.uniprot.org/uniprot/?query=" + \
            str(f) + "+AND+organism:9606&sort=score&limit=1&columns=sequence&format=tab"
        r = requests.get(requestURL, headers={"Accept": "text/plain"})
        if not r.ok:
            r.raise_for_status()
            sys.exit()
        seq1 = r.text.rstrip("\n")
        seq = seq1[9:]
        my_file.write(f"\nLIR Motif Search for {str(set(dict1[f]))}\n")
        LIR_start = None
        for LIRmotif in LIR_motifs:
            if seq.find(LIRmotif) != -1:
                for LIR in re.finditer(LIRmotif, seq):
                    LIR_start = str(LIR.start() + 1)
                    LIR_end = str(LIR.end())
                    my_file.write(
                        LIRmotif +
                        "," +
                        LIR_start +
                        " to " +
                        LIR_end +
                        "\n")
        if not LIR_start:
            my_file.write("No LIR motifs were located in the " +
                          str(set(dict1[f])) + " sequence.")

my_file.close()
