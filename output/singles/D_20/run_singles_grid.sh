#$ -S /bin/bash
#$ -e $JOB_ID.$TASK_ID.e
#$ -o $JOB_ID.$TASK_ID.o
#$ -cwd
#$ -l h_vmem=5G
#$ -l h_cpu=0:30:00
#$ -pe sm 6
#$ -M tomwagg@mpa-garching.mpg.de
#$ -m beas # send an email at begin, ending, abortion and rescheduling of job

export TEMP="/afs/mpa/temp/tomwagg"
export PROJ_DIR="$TEMP/kavli"
export GRID_DIR="$PROJ_DIR/output/singles/D_20"

export MESA_CACHES_DIR="$TEMP/mesa_cache"
export JOB_NAME="$1"

export OMP_NUM_THREADS=$NSLOTS
export MESA_DIR=/afs/mpa/temp/tomwagg/MESA/mesa
export INLIST="inlist_project"

# move to the right directory
cd $GRID_DIR

# ensure we have access to SHMESA
source "$PROJ_DIR/shmesa.sh"

# the range of masses over which we are going to iterate
masses=($(seq 2 0.1 6))
M=${masses[$(($SGE_TASK_ID - 1))]}

echo "Running mass $M"
DIRECTORY=$GRID_DIR/"M_"$M
if [ -d "$DIRECTORY" ]; then echo 'skipping'; exit 0; fi

cp -R ../../../template $DIRECTORY 
cd $DIRECTORY

mesa change $INLIST initial_mass $M
mesa change $INLIST 'xa_central_lower_limit(1)' '1d-5'
mesa change $INLIST min_D_mix 20

./rn

rm star
rm -rf photos

cd -
