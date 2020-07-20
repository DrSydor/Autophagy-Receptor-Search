# Autophagy-Receptor-Search
A Python script that takes in protein ID numbers and names and looks up the amino acid sequence and GO annotations to determine if the protein is a possible autophagy receptor protein.

This script takes a two-column, text delimited file, in which the first column is the protein ID number and the second is the protein name, and returns a text file with only those proteins that possess a putative LIR motif and are annotated (via GO annotation) as ubiquitin-binding. Note that this is a very strict definition of an autophagy receptor protein, so it may miss possible autophagy receptor proteins.
