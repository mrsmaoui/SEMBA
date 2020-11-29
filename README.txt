SEMBA: Single-rEsidue Mutational based Binding Affinity

(c) 2014 Mohamed R. Smaoui, Jerome Waldispuhl
School of Computer Science, McGill University, 3630 University Street, Montreal, Canada
jeromew@cs.mcgill.ca
http://amyloid.cs.mcgill.ca

SEMBA is a program for analyzing the binding affinity of amyloid proteins. 
It uses an energy function that computes the Lennard-Jones, Coulomb, and solvation energies
to determine the effect of a mutation on a protein's stability.

SEMBA runs as command-line application.

INSTALLATION:

SEMBA utilizes the output of several programs to compute stability. In particular, it requires
the prior installation of GROMACS, SCWRL4.

If you have these programs installed, please configure their paths in the file AMLAG.py

If you don't have these programs installed, you will need to download and install them from the
following sources:
GROMACS: http://www.gromacs.org/Downloads
SCWRL4:	http://dunbrack.fccc.edu/scwrl4/ 
PDB2PQR: http://www.ics.uci.edu/~dock/pdb2pqr/userguide.html

Before you use SEMBA, you need to run python install.py to install the AquaSol program that
calculate Solvation energy and the PDB2PQR program. To do that execute the following 2 commands
with the sudo command if needed

	$ sudo tar zxvf AquaSol_Complexes.tgz
	$ sudo ./AquaSol_Complexes/Makefile



USAGE 1: Calculate the binding affinity of a single residue

	$ python SEMBA.py -f protein.pdb -p position -s sequence
		-f protein.pdb: the PDB file you wish to use
		-p position: the residue (amino acid) position you want to calculat the binding
		   affinity for.
		-s sequence: the complete amino acid sequence of the residues in the protein file

	This commands calculates the complete binding affinity of a single residue by computing
	the residue's binding affinity with respect to all 20 amino acids. Each amino acid is 
	paired once on top of the residue being tested. This command returns a table of 20
	runs and calculates the average binding affinity of the residue across all amino acids.


USAGE 2: Calculate the binding affinity of a single residue with resect to a sing amino acid

	$ python SEMBA.py -f protein.pdb -p position -a aminoacid -s sequence
		-f protein.pdb: the PDB file you wish to use
		-p position: the residue (amino acid) position you want to calculat the binding
                   affinity for.
		-a aminoacid: the amino acid letter (ex. K) you wish to calculate the binding
		   affinity for with respect to the residue at position p.
		-s sequence: the complete amino acid sequence of the residues in the protein file

	This commands calculates the binding affinity of a single residue with respect to the
	chosen amino acid. The binding affinity is calculated by analyzing the energy of the
	bonds between the residue and the amino acid placed on top of it.


USAGE 3: Calculate the complete binding affinity of a protein

	$ python SEMBA.py -f protein.pdb -p all -s sequence
		-f protein.pdb: the PDB file you wish to use
		-p position: the residue (amino acid) position you want to calculat the binding
                   affinity for.
		-s sequence: the complete amino acid sequence of the residues in the protein file

	This commands executes 20*n runs where n is the number of residues in the protein. This 
	call explores the binding affinity of all the residues in your chosen protein. This 
	command returns a table with the binding affinity of each residue with respect to each
	of the 20 amino acids, and a table summarizing the average binding affinity of each
	residue.	


