#$ -S /bin/bash
#$ -e $JOB_ID.$TASK_ID.e
#$ -o $JOB_ID.$TASK_ID.o
#$ -cwd
#$ -l h_vmem=5G
#$ -l h_cpu=04:00:00
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
# export GYRE_DIR="$MESA_DIR/gyre/gyre"
export GYRE_DIR="/afs/mpa/temp/tomwagg/kavli/gyre_dev"
export INLIST="inlist_project"

# move to the right directory
cd $GRID_DIR

# the range of masses over which we are going to iterate
masses=($(seq 2 0.1 6))
M=${masses[$(($SGE_TASK_ID - 1))]}

echo "Running mass $M"
DIRECTORY=$GRID_DIR/"M_"$M
cd $DIRECTORY/LOGS

# delete old run
rm -rf *-freqs.dat

find . -name "*.GYRE" -print0 | xargs -0 -P 6 -I{} $PROJ_DIR/GYRE_submitter.sh -i {} -t 1

cd -
