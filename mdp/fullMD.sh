#!/bin/sh

/scratch/msmaoui/gromacs/bin/pdb2gmx -ignh -ff amber99sb-ildn -f $1.pdb -o $1.gro -p $1.top -water tip3p
/scratch/msmaoui/gromacs/bin/editconf -f $1.gro -o $1-PBC.gro -bt triclinic -d 3.0
/scratch/msmaoui/gromacs/bin/grompp -f /scratch/msmaoui/SEMBA/mdp/em-vac-pme.mdp -c $1-PBC.gro -p $1.top -o em-vac.tpr
/scratch/msmaoui/gromacs/bin/mdrun -v -deffnm em-vac

/scratch/msmaoui/gromacs/bin/genbox -cp em-vac.gro -cs spc216.gro -p $1.top -o $1-b4ion.gro
/scratch/msmaoui/gromacs/bin/grompp -f /scratch/msmaoui/SEMBA/mdp/em-sol-pme.mdp -c $1-b4ion.gro -p $1.top -o ion.tpr
echo 13 | /scratch/msmaoui/gromacs/bin/genion -s ion.tpr -o $1-b4em.gro -neutral -conc 0.15 -p $1.top -g ion.log

/scratch/msmaoui/gromacs/bin/grompp -f /scratch/msmaoui/SEMBA/mdp/em-sol-pme.mdp -c $1-b4em.gro -p $1.top -o em-sol.tpr
/scratch/msmaoui/gromacs/bin/mdrun -v -deffnm em-sol

/scratch/msmaoui/gromacs/bin/grompp -f /scratch/msmaoui/SEMBA/mdp/nvt-pr-md.mdp -c em-sol.gro -p $1.top -o nvt-pr.tpr
/scratch/msmaoui/gromacs/bin/mdrun -deffnm nvt-pr

/scratch/msmaoui/gromacs/bin/grompp -f /scratch/msmaoui/SEMBA/mdp/npt-pr-md.mdp -c em-sol.gro -p $1.top -o npt-pr.tpr
/scratch/msmaoui/gromacs/bin/mdrun -deffnm npt-pr

/scratch/msmaoui/gromacs/bin/grompp -f /scratch/msmaoui/SEMBA/mdp/npt-nopr-md.mdp -c em-sol.gro -p $1.top -o npt-nopr.tpr
/scratch/msmaoui/gromacs/bin/mdrun -deffnm npt-nopr

echo 2 2 | /scratch/msmaoui/gromacs/bin/g_rms -s npt-nopr.tpr -f npt-nopr.trr -o $1-bkbone-rmsd.xvg

echo 1 | /scratch/msmaoui/gromacs/bin/trjconv -s npt-nopr.tpr -f npt-nopr.trr -o $1-movie.pdb

# To obtain initial energies, call this
# /scratch/msmaoui/gromacs/bin/g_energy -f npt-nopr.edr -s npt-nopr.tpr -o energies.xvg -b 0 -e 0

#gnuplot> plot "./scwrl_random_2KB8-bkbone-rmsd.xvg" using 1:2 with lines

#Making video in pymol
#dss 
#show cartoon 
#viewport 640,800 
#set ray_trace_frames,1 
#mpng frame_.png
