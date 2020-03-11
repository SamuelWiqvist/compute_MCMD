# Set up for run:

# need this since I use a LU project
#SBATCH -A lu2018-2-22
#SBATCH -p lu

# use gpu nodes
#SBATCH -N 1
#SBATCH -n 1
# #SBATCH --mem-per-cpu=11000
# #SBATCH -C mem256GB

# time consumption HH:MM:SS
#SBATCH -t 150:00:00

# name for script
#SBATCH -J amber_run

# controll job outputs
#SBATCH -o lunarc_output/outputs_amber_run_%j.out
#SBATCH -e lunarc_output/errors_amber_run_%j.err

# notification
#SBATCH --mail-user=samuel.wiqvist@matstat.lu.se
#SBATCH --mail-type=ALL

# load modules
ml load  icc/2017.4.196-GCC-6.4.0-2.28
ml load OpenMPI/2.1.1
ml load Ambertools/18.10

$AMBERHOME/bin/sander -O -i 01_Min.in -o 01_Min.out -p parm7 -c rst7 -r 01_Min.ncrst \ -inf 01_Min.mdinfo
