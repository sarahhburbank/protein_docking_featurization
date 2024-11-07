This repository explores the relatinship between flexbility and performance on a rigid protein-protein docking model using the ABAG, DIPS, and DP5 docking benchmarks. This project shows rhat there is a learnable relationship between flexibility  and performance with a rigid model. This reaffirms that a use case for the rigid heuristic does exist but that ultimately flexible models are required for many docking use cases. 

# benchmark_featurizations folder:

features.py featurizes the proteins based on chemical stucture, size, and evolutionary origin features. There are 3 types of features:

Size Features:
Number of atoms in receptor , Number of atoms in ligand, number of atoms ratio,  Number of unique chains in receptor, Number of unique chains in ligan,number of chains ratio, std dev of residue sizes in receptor,  std dev of residue sizes in ligand

Origin features:
Kingdom of protein origin (eukaryota, viruses, synthetic, bacteria, other, unclassified)
Synthetic Origin Boolean

Chemical Features:
Number of Disulfide Bonds
Number of Predicted Hinge Sites predicted by Packman Hinge 2.8

The results are in benchmark_featurizations.
The order of features in these files are:      

PDB ID, Equidock c-rmsd, hinge_count, kingdom1, synthetic, num_atoms1, num_atoms2, at_ratio,  unique_chains1, unique_chains2,  chain_ratio, std_dev1, std_dev2, disufide_bonds

# ml_analysis folder
This folder uses a random forest classifier to learn the relationship between flexibility and docking performance. The result of this project was to show that there is a strongly learnable relationship between the two.

# benchmark_equidock_crmsd
This folder contains the results of running Equidock, a graph based protein-protein docking model, on the ABAG, DIPS, and DB5 datasets.