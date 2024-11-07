import Bio
from Bio.PDB import PDBParser, Selection, NeighborSearch
import subprocess
from Bio.PDB.PDBParser import PDBParser
from Bio import Entrez
from Bio.SeqUtils import molecular_weight
import os
import numpy as np
import pandas as pd
import urllib
# Create a PDBParser object
parser = PDBParser()
def size_features(file_name):
    #try: 
    structure = parser.get_structure(file_name, f'/Users/sarahburbank/Desktop/iw/pdb_analysis/ftr/proteins/{file_name}')
    
    #----------------size variables--------------------#
    #num atoms
    all_atoms = structure.get_atoms()
    #atomic_weight = 0
    # Get the atom type

    num_atoms = sum(1 for __ in all_atoms)
    #print("atom count: ", num_atoms)


    print( num_atoms)
size_features('3szk.pdb')