# Set up for run:

# need this since I use a LU project
#SBATCH -A snic2019-3-630


# use gpu nodes
#SBATCH -N 1
#SBATCH -n 1
# #SBATCH --mem-per-cpu=11000
# #SBATCH -C mem256GB

# time consumption HH:MM:SS
#SBATCH -t 1:00:00

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

mpirun -np 16 $AMBERHOME/bin/sander -O -i heat1.in -p TC5b.prmtop -c min1.ncrst -r heat1.ncrst -o heat1.out -x heat1.nc
gzip -9 heat1.nc
mpirun -np 16 $AMBERHOME/bin/sander -O -i heat2.in -p TC5b.prmtop -c heat1.ncrst -r heat2.ncrst -o heat2.out -x heat2.nc
gzip -9 heat2.nc
mpirun -np 16 $AMBERHOME/bin/sander -O -i heat3.in -p TC5b.prmtop -c heat2.ncrst -r heat3.ncrst -o heat3.out -x heat3.nc
gzip -9 heat3.nc
mpirun -np 16 $AMBERHOME/bin/sander -O -i heat4.in -p TC5b.prmtop -c heat3.ncrst -r heat4.ncrst -o heat4.out -x heat4.nc
gzip -9 heat4.nc
mpirun -np 16 $AMBERHOME/bin/sander -O -i heat5.in -p TC5b.prmtop -c heat4.ncrst -r heat5.ncrst -o heat5.out -x heat5.nc
gzip -9 heat5.nc
mpirun -np 16 $AMBERHOME/bin/sander -O -i heat6.in -p TC5b.prmtop -c heat5.ncrst -r heat6.ncrst -o heat6.out -x heat6.nc
gzip -9 heat6.nc
mpirun -np 16 $AMBERHOME/bin/sander -O -i heat7.in -p TC5b.prmtop -c heat6.ncrst -r heat7.ncrst -o heat7.out -x heat7.nc
gzip -9 heat7.nc
