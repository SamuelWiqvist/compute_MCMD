#!/bin/sh

# Set up for run:

# need this since I use a LU project
#SBATCH -A lu2019-2-19
#SBATCH -p lu

# use gpu nodes
#SBATCH -N 1
#SBATCH -n 1


# time consumption HH:MM:SS
#SBATCH -t 1:00:00

# name for script
#SBATCH -J amber_run

# controll job outputs
#SBATCH -o lunarc_outputs_amber_md%j.out
#SBATCH -e lunarc_outputs_amber_md%j.err

# notification
#SBATCH --mail-user=samuel.wiqvist@matstat.lu.se
#SBATCH --mail-type=ALL

# load modules
ml load  icc/2017.4.196-GCC-6.4.0-2.28
ml load OpenMPI/2.1.1
ml load Ambertools/18.10

# Run minimization
$AMBERHOME/bin/sander -O -i 01_Min.in -o 01_Min.out -p parm7 -c rst7 -r 01_Min.ncrst \
-inf 01_Min.mdinfo

# Run heating MD
$AMBERHOME/bin/sander -O -i 02_Heat.in -o 02_Heat.out -p parm7 -c 01_Min.ncrst \
-r 02_Heat.ncrst -x 02_Heat.nc -inf 02_Heat.mdinfo

# Run production MD
$AMBERHOME/bin/sander -O -i 03_Prod.in -o 03_Prod.out -p parm7 -c 02_Heat.ncrst \
-r 03_Prod.ncrst -x 03_Prod.nc -inf 03_Prod.info 
