#!/bin/bash
source ../../shmesa.sh

INLIST="inlist_project"
OUTPUT="./"

single () {
    local M=$1
    local MT=$2
    
    PARAMS="MT_"$MT
    echo "Running mass $M" with $PARAMS
    DIRECTORY=$OUTPUT/"M_"$M"-"$PARAMS
    if [ -d "$DIRECTORY" ]; then echo 'skipping'; return 0; fi
    
    cp -R ../../template $DIRECTORY 
    cd $DIRECTORY 
    
    if [[ $MT == 1 ]]; then
        mesa change $INLIST use_other_adjust_mdot ".true."
        mesa change $INLIST 'x_logical_ctrl(1)' ".true."
        mesa change $INLIST 'x_ctrl(1)' 0.65
        mesa change $INLIST 'x_ctrl(2)' 3.5
        mesa change $INLIST 'x_ctrl(3)' 2.5d-7
    fi
    mesa change $INLIST initial_mass $M
    
    ./rn

    rm star
    rm -rf photos
    
    cd -
}

run_gyre () {
    local M=$1
    local os=$2

    PARAMS="os_"$os
    echo "Running GYRE for run M=$M" with $PARAMS
    DIRECTORY=$OUTPUT/"M_"$M"-"$PARAMS
    if [ -d "$DIRECTORY" ]; then echo 'skipping'; return 0; fi

    cd $DIRECTORY/LOGS

    find . -name "*.GYRE" -print0 | xargs -0 -P 6 -I{} ~/Documents/research/kavli-2023/algol-seismology/gyre6freqs.sh -i {} -t 1

    cd -
}

single 3 0
single 4 1
# single 3 1
single 3.51317 0
