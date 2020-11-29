#!/usr/bin/env python
from subprocess import call
import string
import sys
import AquasolCPLX_routineBinding

def calculateE(filename):
	#Calculate Coulomb and LJ energies
	call("cd /scratch/msmaoui/MeasureAffinity/affinity/; /scratch/msmaoui/FibrilMutant/Amylin/mdp/zeroMD.sh %s > /scratch/msmaoui/MeasureAffinity/affinity/trash/result_%s.tmp 2>&1" % (filename,filename), shell=True)

        call("tail -10 /scratch/msmaoui/MeasureAffinity/affinity/trash/result_%s.tmp > /scratch/msmaoui/MeasureAffinity/affinity/trash/resultshort_%s.tmp" % (filename, filename), shell=True)

        fileGrid = open('/scratch/msmaoui/MeasureAffinity/affinity/trash/resultshort_%s.tmp' % filename, 'r')
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
        energy = AquasolCPLX_routineBinding.applyAquaSol(temp, filename+'.pdb', filename)
	(Fw, Fv, Nw) = energy.split( )


	call("rm -f /scratch/msmaoui/MeasureAffinity/affinity/trash/resultshort_%s.tmp /scratch/msmaoui/MeasureAffinity/affinity/trash/result_%s.tmp" % (filename, filename), shell=True)
        call("rm -f /scratch/msmaoui/MeasureAffinity/affinity/%s*" % filename, shell=True)
        call("rm -f /scratch/msmaoui/MeasureAffinity/affinity/#*", shell= True)


	return(Coul*0.238845897, LJ*0.238845897, Fw, Fv, Nw)


def mutateSequence(x, m):
	originalSeq = list("KCNTATCATQRLANFLVHSSNNFGAILSSTNVGSNTY")
	originalSeq[x-1] = m
	return "".join(originalSeq)	
		
def removeResidues(x, filename):
	fileInfo = open('/scratch/msmaoui/MeasureAffinity/affinity/%s' % filename, 'r')
        fileLines = fileInfo.readlines()
	
	fileInfo_n = open('/scratch/msmaoui/MeasureAffinity/affinity/tmp_%s' % filename, 'w')
	
	for one_line in fileLines:
		ones = string.split(one_line)
		if(len(ones) > 8):
			if(int(ones[5]) < 371):
				fileInfo_n.write(one_line)
			if(int(ones[5]) == (370+int(x))):
				fileInfo_n.write(one_line)
			
	fileInfo_n.close()
	fileInfo.close()
	call('rm /scratch/msmaoui/MeasureAffinity/affinity/%s' % filename, shell=True)
	call('mv /scratch/msmaoui/MeasureAffinity/affinity/tmp_%s /scratch/msmaoui/MeasureAffinity/affinity/%s' % (filename, filename), shell=True)

def main():
	mutationList = list("ARNDCQEGHILKMFPSTWYV")

	count = 0
	print "Amylin Binding Affinity Landscape"
	print "Count,P,M,Coul1,LJ1,Fv1,Fw1,Nw1,Coul2,LJ2,Fv2,Fw2,Nw2"

	for x in range(1, 38):
		for m in mutationList:
			#Update Counter
			count = count + 1

			if(count > 22):
				#Perform Mutation on original sequence
				sequence = mutateSequence(x, m)
	
				#try:	
				#Generate .fasta file for interaction with Amyloid
				call("echo KCNTATCATQRLANFLVHSSNNFGAILSSTNVGSNTY%s > /scratch/msmaoui/MeasureAffinity/affinity/2kb8_aff_%d-%s.fasta" % (sequence, x, m), shell=True)
				#Generate Structure file by SCWRL for interaction with Amyloid
				call("/scratch/msmaoui/SCWRL4.0/Scwrl4 -i /scratch/msmaoui/MeasureAffinity/2kb8_dimer.pdb -o /scratch/msmaoui/MeasureAffinity/affinity/scwrl_2kb8_aff_%d-%s.pdb -s /scratch/msmaoui/MeasureAffinity/affinity/2kb8_aff_%d-%s.fasta >/dev/null" % (x, m, x, m), shell=True)

				#remove all residues but the mutated one
				removeResidues(x, "scwrl_2kb8_aff_%d-%s.pdb" % (x, m))

				(Coul1, LJ1, Fw1, Fv1, Nw1) = calculateE("scwrl_2kb8_aff_%d-%s" % (x, m))
	
				call("echo KCNTATCATQRLANFLVHSSNNFGPILPPTNVGSNTY%s > /scratch/msmaoui/MeasureAffinity/affinity/Pram_aff_%d-%s.fasta" % (sequence, x, m), shell=True)
        	               	#Generate Structure file by SCWRL for interaction with Amyloid
                	       	call("/scratch/msmaoui/SCWRL4.0/Scwrl4 -i /scratch/msmaoui/MeasureAffinity/2kb8_dimer.pdb -o /scratch/msmaoui/MeasureAffinity/affinity/scwrl_Pram_aff_%d-%s.pdb -s /scratch/msmaoui/MeasureAffinity/affinity/Pram_aff_%d-%s.fasta >/dev/null" % (x, m, x, m), shell=True)

                       		#remove all residues but the mutated one
                       		removeResidues(x, "scwrl_Pram_aff_%d-%s.pdb" % (x, m))

                       		(Coul2, LJ2, Fw2, Fv2, Nw2) = calculateE("scwrl_Pram_aff_%d-%s" % (x, m))


				print str(count) + "," + str(x) + "," + m + "," + str(Coul1) + "," + str(LJ1) + "," + str(Fw1) + "," + str(Fv1) + "," + str(Nw1) + "," + str(Coul2) + "," + str(LJ2) + "," + str(Fw2) + "," + str(Fv2) + "," + str(Nw2)

				#except:
				#	print str(count) + "," + str(x) + "," + m + ",0.0,0.0,0.0,0.0,0.0," + sequence

				sys.stdout.flush()

main()
