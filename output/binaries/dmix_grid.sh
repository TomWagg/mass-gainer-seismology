#$ -S /bin/bash
#$ -e $JOB_ID.$TASK_ID.e
#$ -o $JOB_ID.$TASK_ID.o
#$ -cwd
#$ -l h_vmem=10G
#$ -l h_cpu=3:00:00
#$ -pe sm 6
#$ -M tomwagg@mpa-garching.mpg.de
#$ -m beas # send an email at begin, ending, abortion and rescheduling of job

export TEMP="/afs/mpa/temp/tomwagg"
export PROJ_DIR="$TEMP/kavli"
export GRID_DIR="$PROJ_DIR/output/binaries"

export MESA_CACHES_DIR="$TEMP/mesa_cache"
export JOB_NAME="$1"

export OMP_NUM_THREADS=$NSLOTS
export MESA_DIR=/afs/mpa/temp/tomwagg/MESA/mesa

# move to the right directory
cd $GRID_DIR

# ensure we have access to SHMESA
source "$PROJ_DIR/shmesa.sh"

# the range of masses over which we are going to iterate
dmixs=($(seq 15 5 45))
D=${dmixs[$(($SGE_TASK_ID - 1))]}

echo "Running mass $M"
DIRECTORY=$GRID_DIR/"D_"$D
if [ -d "$DIRECTORY" ]; then echo 'skipping'; exit 0; fi

cp -R ../../template_binary $DIRECTORY 
cd $DIRECTORY

mesa change "inlist1" min_D_mix $D
mesa change "inlist2" min_D_mix $D

./rn

cd -
