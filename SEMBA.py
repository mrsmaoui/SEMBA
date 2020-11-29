#!/usr/bin/env python

############# MODIFY PROGRAM PATHS ACCORDINGLY #############
SEMBA_PATH="/scratch/msmaoui/SEMBA"
SCWRL4_PATH="/scratch/msmaoui/SCWRL4.0"
GROMACE_BIN_PATH="/scratch/msmaoui/gromacs/bin"
PDB2PQR_PATH="/scratch/msmaoui/pdb2pqr-1.7"
AQUASOL_PATH="/data/pkgs/AquaSol_Complexes/bin"
############# YOU DONT NEED TO MODIFY PAST THIS POINT #############


from subprocess import call
import string
import sys 
import AquasolCPLX_routine
import time
import os

def main(argv):
	error = "ERROR: Please consult the README.txt FILE in your SEMBA directory to properly use this program"

	#verifying input variables are correct
	for i in range(len(argv)):
        	if( i % 2 == 1):        #check odd values
                	#exit if command doesn't exist
                	if not (argv[i] == '-f' or argv[i] == '-p' or argv[i] == '-a' or argv[i] == '-s'):
                	        print error
                	        sys.exit(0)

	#get input parameters
	filename = ''
	position = ''
	mutation = ''
	sequence = ''

	try:
        	for i in range(len(argv)):
        	        if(argv[i] == '-f'): filename = argv[i+1]
        	        elif (argv[i] == '-p'): position = argv[i+1]
        	        elif (argv[i] == '-a'): mutation = argv[i+1]
			elif (argv[i] == '-s'): sequence = argv[i+1]

	except Exception:
	        print error
        	sys.exit(0)

	#check that input file is of pdb format
	if not filename[len(filename)-4:len(filename)] == '.pdb':
	        print error
		print "File entered does not have a .pdb extension"
        	sys.exit(1)

	filename = filename[0:(len(filename)-4)]

	if(sequence == ""):
		print error
		print "You didn't specify the amino acid sequence with parameter -s"
		sys.exit(1)

	if(mutation != "" and len(mutation) != 1):
		mutation = mutation.upper()
		#check that the mutation is a valid character
		mutationList = list("ARNDCQEGHILKMFPSTWYV")
		if(mutation not in mutationList):
			print "ERROR: Incorrect mutation choice"
			sys.exit(0)	
	

	#Determine which USAGE the user intends to perform
	if( (position != "all") and mutation == ''):
		params = {"USAGE": 1, "filename":filename, "position":position, "mutation":'', "sequence":sequence}
		return USAGE1(params)

	if( (position != 'all') and mutation != ""):
		params = {"USAGE":2, "filename":filename, "position":position, "mutation":mutation, "sequence":sequence}
		return USAGE2(params)

	if( (position == 'all') ):
                params = {"USAGE":3, "filename":filename, "position":position, "mutation":'', "sequence":sequence}        
                return USAGE3(params)

	else:
		print error
		sys.exit(0)
	



def calculateE(usage_path, filename):
	#Calculate Coulomb and LJ energies
	call("cd %s; %s/mdp/zeroMD.sh %s %s %s > %s/trash/result_%s.tmp 2>&1" % (usage_path, SEMBA_PATH, GROMACE_BIN_PATH, filename, SEMBA_PATH, usage_path, filename), shell=True)

        call("tail -10 %s/trash/result_%s.tmp > %s/trash/resultshort_%s.tmp" % (usage_path, filename, usage_path, filename), shell=True)

        fileGrid = open('%s/trash/resultshort_%s.tmp' % (usage_path, filename), 'r')
        fileLines = fileGrid.readlines()

        LJ = 0.0
        Coul = 0.0
        for one_line in fileLines:
	        ones = string.split(one_line)
                if(len(ones) > 1 and ones[0] == 'LJ-14'):
                        LJ = float(ones[1])
                if(len(ones) > 1 and ones[0] == 'Coulomb'):
                        Coul = float(ones[2])
	
	#Calculate solvation energy
	temp = 300
        energy = AquasolCPLX_routine.applyAquaSol(AQUASOL_PATH, usage_path, PDB2PQR_PATH, SEMBA_PATH, temp, filename+'.pdb', filename)
	(Fw, Fv, Nw) = energy.split( )


	call("rm -f %s/trash/resultshort_%s.tmp %s/trash/result_%s.tmp" % (usage_path, filename, usage_path, filename), shell=True)


	return(Coul*0.238845897, LJ*0.238845897, Fw, Fv, Nw)


def mutateSequence(x, m, sequence):
	originalSeq = list(sequence)
	originalSeq[x-1] = m
	return "".join(originalSeq)	
	
def removeResidues(x, path, filename, seq_length):

        fileInfo = open('%s/%s' % (path, filename), 'r')
        fileLines = fileInfo.readlines()

        fileInfo_n = open('%s/tmp_%s' % (path, filename), 'w')

	starting_res = -1
        for one_line in fileLines:
                ones = string.split(one_line)
                if(len(ones) > 8):
			if(starting_res == -1):
				starting_res = int(ones[5])

                        if(int(ones[5]) < (starting_res + seq_length)):
                                fileInfo_n.write(one_line)
                        if(int(ones[5]) == (starting_res + seq_length + int(x))):
                                fileInfo_n.write(one_line)

        fileInfo_n.close()
        fileInfo.close()
        call('rm %s/%s' % (path, filename), shell=True)
        call('mv %s/tmp_%s %s/%s' % (path, filename, path, filename), shell=True)



def USAGE1(params):

	filename = params['filename']
	position = int(params['position'])
	sequence = params['sequence']

	mutationList = list("ARNDCQEGHILKMFPSTWYV")

	#Calculate the energy of a protein before binding
	beforebinding = simulate('', sequence, -1, filename)

	#Create a second protein on top of the amyloid protein: dimer protein
	filename2 = createDimer(filename)

	for mutation in mutationList:
		#Perform Mutation on original sequence
		mutant_seq = mutateSequence(position, mutation, sequence)
		afterbinding = simulate(mutation, mutant_seq, position, filename2)


	


def USAGE2(params):

	filename = params['filename']
        mutation = params['mutation'].split(',')
        position = [int(x) for x in params['position'].split(',')]
        sequence = params['sequence']

	mutant_seq = sequence
	for m in range(0,len(position)):
		mutant_seq = mutateSequence(position[m], mutation[m], mutant_seq)

	print "Position,Mutation,Coul,LJ,Solvation,Sequence"
        print simulate('multp', mutant_seq, 0, filename)

        mutation = "none"
        print simulate(mutation, sequence, 0, filename)


def USAGE3(params):
	filename = params['filename']
        position = int(params['position'])
        sequence = params['sequence']

	mutationList = list("ARNDCQEGHILKMFPSTWYV")
	
	print "Position,Mutation,Coul,LJ,Solvation,Sequence"
	for mutation in mutationList:
		mutant_seq = mutateSequence(position, mutation, sequence)
		print simulate(mutation, mutant_seq, position, filename)


def simulate(mutation, sequence, position, filename):
	os.chdir(SEMBA_PATH)	
	timestr = time.strftime("%Y%m%d-%H%M%S")
	usage_path = SEMBA_PATH + '/tmp/' + timestr

	#Create a workspace directory
	call('mkdir %s' % usage_path, shell=True)
	call('mkdir %s/trash' % usage_path, shell=True)

	filename2 = createDimer(filename)

	#Generate .fasta file 
	call('echo %s > %s/%s_%d-%s.fasta' % (sequence, usage_path, filename,position,mutation), shell=True)
	#Perform mutation and generate structure file
	call("%s/Scwrl4 -i %s/%s.pdb -o %s/scwrl_%s_%d-%s.pdb -s %s/%s_%d-%s.fasta >/dev/null" % (SCWRL4_PATH, SEMBA_PATH, filename, usage_path, filename, position, mutation, usage_path, filename, position, mutation), shell=True)	

	removeResidues(position, usage_path, "scwrl_%s_%d-%s.pdb" % (filename, position, mutation), len(sequence))

	(Coul, LJ, Fw, Fv, Nw) = calculateE(usage_path, "scwrl_%s_%d-%s" % (filename, position, mutation))

	solvation = float(Fw) - float(Fv) - (float(Nw) * -0.0352669760297406)

	#remove directory
	call('rm -rf %s/tmp/%s' % (SEMBA_PATH, timestr), shell=True)
	
	return str(position) + "," + mutation + "," + str(Coul) + "," + str(LJ) + "," + str(solvation) + "," + sequence


	
def USAGE4(params):
	filename = params['filename']
	sequence = params['sequence']
	(start_res, end_res) = params['position'].split(',')
	start_res = int(start_res)
	end_res = int(end_res)

	mutationList = list("ARNDCQEGHILKMFPSTWYV")

	count = 0
	print "Position,Mutation,Coul,LJ,Solvation,Sequence"

	
	for x in range(start_res, (end_res+1)):
		for m in mutationList:
			#Update Counter
			count = count + 1

			#Perform Mutation on original sequence
			mutant_seq = mutateSequence(x, m, sequence)
	
			print simulate(m, mutant_seq, x, filename)
	
			sys.stdout.flush()
main(sys.argv)

