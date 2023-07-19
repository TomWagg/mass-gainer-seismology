#$ -S /bin/bash
#$ -e $JOB_ID.e
#$ -o $JOB_ID.o
#$ -cwd
#$ -l h_vmem=10G
#$ -l h_cpu=06:00:00
#$ -l qname=hilbert
#$ -pe sm 12
#$ -M tomwagg@mpa-garching.mpg.de
#$ -m beas # send an email at begin, ending, abortion and rescheduling of job

export TEMP="/afs/mpa/temp/tomwagg"
export PROJ_DIR="$TEMP/kavli"
export MESA_CACHES_DIR="$TEMP/mesa_cache"
export JOB_NAME="$1"

export OMP_NUM_THREADS=$NSLOTS
export MESA_DIR=/afs/mpa/temp/tomwagg/MESA/mesa

./rn












