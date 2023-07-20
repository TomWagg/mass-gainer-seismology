#$ -S /bin/bash
#$ -e $JOB_ID.$TASK_ID.e
#$ -o $JOB_ID.$TASK_ID.o
#$ -cwd
#$ -l h_vmem=10G
#$ -l h_cpu=0:30:00
#$ -pe sm 10
#$ -M tomwagg@mpa-garching.mpg.de
#$ -m beas # send an email at begin, ending, abortion and rescheduling of job

export TEMP="/afs/mpa/temp/tomwagg"
export PROJ_DIR="$TEMP/kavli"
export GRID_DIR="$PROJ_DIR/output/cluster/singles_grid"

export MESA_CACHES_DIR="$TEMP/mesa_cache"
export JOB_NAME="$1"

export OMP_NUM_THREADS=$NSLOTS
export MESA_DIR=/afs/mpa/temp/tomwagg/MESA/mesa

# move to the right directory
cd $GRID_DIR

# ensure we have access to SHMESA
source "$PROJ_DIR/shmesa.sh"

# the range of masses over which we are going to iterate
masses=($(seq 3 0.1 6))
TASK_ID=2 # TODO REMOVE THIS
M=${masses[$TASK_ID]}

echo "Running mass $M"
DIRECTORY=$GRID_DIR/"M_"$M
if [ -d "$DIRECTORY" ]; then echo 'skipping'; return 0; fi

cp -R ../../../template $DIRECTORY 
cd $DIRECTORY 

mesa change $INLIST initial_mass $M

./rn

rm star
rm -rf photos

cd -
