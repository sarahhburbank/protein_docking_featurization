from Bio.PDB import PDBParser, NeighborSearch
import subprocess
from Bio.PDB.PDBParser import PDBParser
from Bio import Entrez
from Bio.SeqUtils import molecular_weight
import os
import numpy as np
import pandas as pd
import urllib
import csv

parser = PDBParser()
Entrez.email = "sb9342@princeton.edu"

#can return eukaryota, viruses, synthetic, bacteria, other, unclassified"
def origin_features(id):
    #get correct type of pdb#
    #new 
    url = f'https://files.rcsb.org/download/{id}.pdb'
    with urllib.request.urlopen(url) as response:
        pdb_data = response.read().decode('utf-8')    
    try:
        lines = pdb_data.split('\n')
        tax_ids = set()
        kingdoms = set()
        synthetic = False
        enzyme = False
        for line in lines:
            # Split the line into fields using whitespace as delimiter
            fields = line.split()
            if len(fields) >= 4 and "SYNTHETIC" in fields[2]:
                synthetic = True
            if len(fields) >= 4 and "ENZYME" in fields[2]:
                enzyme = True
            # Check if there are at least three fields and if the third field contains "ORGANISM_TAXID"
            if len(fields) >= 4 and "ORGANISM_TAXID" in fields[2]:
                # Append the line to the list of organism taxid rows
                tax_ids.add(fields[3]) 
                #print("chain taxid", line)
                handle = Entrez.efetch(db="taxonomy", id=str(fields[3]), retmode="xml")
                records = Entrez.read(handle)
                handle.close()
                for record in records:
                    for line in record['LineageEx']:
                        if line['Rank'] == 'superkingdom' and line['ScientificName'] not in kingdoms:
                            kingdoms.add(line['ScientificName'])
                            #print("sci name",line['ScientificName']) 
        if len(kingdoms)  == 0:
            kingdoms.add('other') 
    except:
        return 0,0        
    return kingdoms, synthetic, 



# Parse the PDB file
def size_features(file_name):
    structure = parser.get_structure(file_name, f'/Users/sarahburbank/Desktop/iw/abag/{file_name}')
    
    #----------------size variables--------------------#
    #num atoms
    all_atoms = structure.get_atoms()
    #atomic_weight = 0
    # Get the atom type

    num_atoms = sum(1 for __ in all_atoms)
    #print("atom count: ", num_atoms)

    #num stuff
    residue_counts = {}
    unique_chains = 0
    total_chains = 0
    chain_ids = set()
    residue_sizes = []

    for chain in structure[0]:
        #print("chain", chain)
        #----------------size information--------------------#
        total_chains+=1

        if chain.id not in chain_ids:
            chain_ids.add(chain.id)
        res_chain = 0
        tot_residues = 0
        for residue in chain:
            #print("res", residue)
            res = residue.get_resname()
            #atomic_weight += molecular_weight(res, "protein")
            res_chain+=1
            tot_residues+= 1
            '''residue_name = residue.get_resname()
            if residue_name in residue_counts:
                residue_counts[residue_name] += 1
            else:
                residue_counts[residue_name] = 1'''
        residue_sizes.append(tot_residues)
        unique_chains = len(chain_ids)
        std_dev_res_chain = np.std(residue_sizes)
    return num_atoms, residue_counts, unique_chains, total_chains, std_dev_res_chain 
    
    # Access information about the structure
def chemical(name):
    structure = parser.get_structure(name, f'/Users/sarahburbank/Desktop/iw/pdb_analysis/ftr/pdb_files/{name}')
    cysteines = [residue for residue in structure.get_residues() if residue.get_resname() == "CYS"]
    disulfide_bonds = 0

    atom_list = [atom for residue in cysteines for atom in residue if atom.element == "S"]
    searcher = NeighborSearch(atom_list)
    for atom in atom_list:
        close_atoms = searcher.search(atom.coord, 2.2, "R")  # "R" means search radius
        for close_atom in close_atoms:
            if close_atom != atom:  # ensure it's not the atom itself
                pair = tuple(sorted([atom, close_atom], key=lambda x: x.serial_number))
                if pair not in disulfide_bonds:
                    disulfide_bonds += 1
    return disulfide_bonds

def two_proteins(filename_1, filename_2, id ):
    #-----------------------size features------------------------#
    num_atoms1, residue_counts1, unique_chains1, total_chains1, std_dev1 = size_features(filename_1)
    num_atoms2, residue_counts2, unique_chains2, total_chains2, std_dev2 = size_features(filename_2)
        
    ratio1 = num_atoms1 / num_atoms2
    ratio2 = num_atoms2 /num_atoms1
    if ratio1 > ratio2:
        at_ratio = ratio1
    else:
        at_ratio = ratio2
    
    chain_ratio1 = total_chains1 / total_chains2
    chain_ratio2 = total_chains2 / total_chains1
    if chain_ratio1 > chain_ratio2:
        chain_ratio = chain_ratio1
    else:
        chain_ratio = chain_ratio2
    #-----------------------origin features------------------------#
    kingdoms1, synth1 = origin_features(id)
    
    if synth1 == True:
        synth1_num = 1
    else:
        synth1_num = 0
    
    if type(kingdoms1) != int :
       
        for thing in kingdoms1:
            if thing == 'Eukaryota':
                king1 = (1, 0 , 0, 0)
            elif thing == 'Viruses':
                king1 = (0, 1, 0, 0)
            elif thing == 'Bacteria':
                king1 = (0, 0, 1, 0)
            elif thing == 'Other' or thing == 'Unclassified':
                king1 = (0, 0, 0, 1)
            else: 
                king1 = (0,0,0,0) 
                #print('king1', king1) 
            '''if len(king1 > 1):
                king1 = (king1[0])'''
    else:
            king1 = (0,0,0,0) 
            print("returned 0,0")
            
    
    #----------------------- feature row ------------------------#
    #can include residue counts later if necessary
    #feature order:
    '''name, confidence score/label	origin1, origin2,	synth 1 	synth 2 	
    #atoms 1 	#atoms 2 	atom num ratio	atom weight 1 	atom weight 2	atomic weight ratio	
    # unique chains1	unique chains 2	total chains 1 	total chains 2, chain ratio
    # std dev res/chain 1	std dev res/chain 2 	 '''
    #atomic_weight1, atomic_weight2, at_ratio_w
    features =  king1, synth1_num, num_atoms1, num_atoms2, at_ratio,  unique_chains1, unique_chains2, total_chains1, total_chains2, chain_ratio, std_dev1, std_dev2
    return features

import csv
def make_tensor_single():
    for pdbfile in os.listdir('/Users/sarahburbank/Desktop/iw/pdb_analysis/ftr/pdb_files'):
        print(pdbfile)
        num_atoms, residue_counts, unique_chains, total_chains, std_dev_res_chain,  = size_features(pdbfile)
        kingdoms, synthetic,  = origin_features(pdbfile)
        row = [pdbfile, kingdoms,  synthetic, num_atoms, residue_counts, unique_chains, total_chains, std_dev_res_chain ]
        with open('/Users/sarahburbank/Desktop/iw/pdb_analysis/ftr/result.csv', 'a', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(row)

import re
def hinge(ligand_file_name):
    if ligand_file_name[2] == '_':#dips 
            parts = ligand_file_name.split("_")  # Split the string by "_"
            id = parts[1].split(".")[0] #gives 4aqn
            pdb_file_hinge_parts =  ligand_file_name.split(".")
            pdb_file_name = pdb_file_hinge_parts[1] #gives aq_4aqn

    elif ligand_file_name[4]== '_':#db5
            parts = ligand_file_name.split("_")  # Split the string by "_"
            id = parts[0]
            pdb_file_name = id #might be wrong 
    
    try:
        print("trying  hinge")
        
        command = f'packman hinge 2.8 {pdb_file_name}.pdb --pdbid {id} --e_clusters 4 --minhnglen 5'
        print(command)
        result = subprocess.run(command, shell=True, text=True, capture_output=True)
        hinge_count = result.stdout.count('Hinge #')
        print("hinge count", hinge_count)
        return hinge_count
       
            
    except:
        print('-1')
        return -1
    
def many_hinge():
    #code input :packman hinge 2.8 aq_4aqn.pdb1_0.dill_l_b.pdb --pdbid 4aqn --e_clusters 4 --minhnglen 5
    #hing count 0
    #hand written input packman hinge 2.8 aq_4aqn.pdb --pdbid 4aqn --e_clusters 4 --minhnglen 5
    all_crmsd_df = pd.read_csv('/Users/sarahburbank/Desktop/iw/pdb_analysis/ftr/all_crmsd.csv',header=0, names=['ligand', 'receptor', 'crmsd'])

    for row in all_crmsd_df.itertuples():
        ligand_file_name = row[1]
        # example dips aq_4aqn.pdb1_0.dill_l_b.pdb
        if ligand_file_name[2] == '_':#dips 
            parts = ligand_file_name.split("_") 
            id = parts[1].split(".")[0] 
            pdb_file_hinge_parts =  ligand_file_name.split(".")
            pdb_file_name = pdb_file_hinge_parts[1] #gives aq_4aqn

        elif ligand_file_name[4]== '_':#db5
            parts = ligand_file_name.split("_")
            id = parts[0]
            pdb_file_name = id 
            
        hinge(id, pdb_file_name)
         
def make_two_protein_tensor(): 
    all_crmsd_df = pd.read_csv('/Users/sarahburbank/Desktop/iw/pdb_analysis/ftr/abag_crmsd.csv',header=0, names=['ligand', 'receptor', 'crmsd','irmsd'])
    
    #for row in all_crmsd_df.itertuples():
    for index, row in enumerate(all_crmsd_df.itertuples(), start=1):
        if index < 84:
            continue
        ligand_file_name = row[1]
        receptor_file_name = row[2]
        print("receptor name",receptor_file_name)
        crmsd = row[3]
        
        if ligand_file_name[2] == '_':#dips 
            parts = ligand_file_name.split("_")  # Split the string by "_"
            id = parts[1].split(".")[0]
            print("id", id)

        elif ligand_file_name[4]== '_':#db5
            parts = ligand_file_name.split("_")  # Split the string by "_"
            id = parts[0]
            print("id", id)
        else: print(f"FILE NAME ERROR FOR {ligand_file_name}")
        hinge_count= hinge(ligand_file_name)
        king1, synth1_num, num_atoms1, num_atoms2, at_ratio,  unique_chains1, unique_chains2, total_chains1, total_chains2,chain_ratio, std_dev1, std_dev2 = two_proteins(ligand_file_name, receptor_file_name,id )
        row = [id, crmsd, hinge_count, king1, synth1_num, num_atoms1, num_atoms2, at_ratio,  unique_chains1, unique_chains2,  chain_ratio, std_dev1, std_dev2]
        with open('/Users/sarahburbank/Desktop/iw/pdb_analysis/ftr/abag_features.csv', 'a', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(row)



def calculate_disulfide_ratio(pdb_filename):
    parser = PDBParser()
    structure = parser.get_structure('structure', pdb_filename)

    disulfide_bonds = 0
    total_atoms = 0

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.resname == 'CYS' and residue.id[0] == ' ':
                    if 'SG' in residue:
                        if residue['SG'].get_parent().id[0] == ' ':
                            disulfide_bonds += 1
                    else: return 0
                total_atoms += len(residue)

    return disulfide_bonds / total_atoms if total_atoms != 0 else 0

csv_filename = '/Users/sarahburbank/Desktop/iw/pdb_analysis/ftr/ml_analysis/all_features .csv'
output_csv_filename = '/Users/sarahburbank/Desktop/iw/pdb_analysis/ftr/features_disulfide.csv'

with open(csv_filename, 'r') as csvfile,  open(output_csv_filename, 'w', newline='') as outfile:
    csvreader = csv.reader(csvfile)
    csvwriter = csv.writer(outfile)

    for index, row in enumerate(csvreader):
        ligand = row[0]
        id = ligand[:4]
        print(id)
        if os.path.exists(f'/Users/sarahburbank/Desktop/iw/all_proteins/{id}_r_u.pdb'):
            pdb = f'/Users/sarahburbank/Desktop/iw/all_proteins/{id}_r_u.pdb'
        else:
            pdb = f'/Users/sarahburbank/Desktop/iw/all_proteins/{id}_l_u.pdb'
        ratio1 = calculate_disulfide_ratio(pdb)
        print(ratio1)
        row1 = [f'{id}', ratio1]
        csvwriter.writerow(row1)
