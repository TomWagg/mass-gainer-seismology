#!/bin/bash

#### Converter for .GYRE file to oscillation mode frequencies with GYRE 
#### Author: Earl Bellinger ( bellinger@phys.au.dk ) 
#### Max Planck Institute for Astrophysics, Garching, Germany 
#### Stellar Astrophysics Centre, Aarhus University, Denmark 

#### (Editted by Tom Wagg for very specific project things)

### Parse command line tokens 
HELP=0
EIGENF=0
SAVE=0
RADIAL=0
FGONG=0
OMP_NUM_THREADS=1
LOWER=0.25
UPPER=10
CONVERT=0
UNITS="CYC_PER_DAY"
RESOLUTION=0
DIPOLE=0

while [ "$#" -gt 0 ]; do
  case "$1" in
    -h) HELP=1; break;;
    -i) INPUT="$2"; shift 2;;
    -o) OUTPUT="$2"; shift 2;;
    -t) OMP_NUM_THREADS="$2"; shift 2;;
    -l) LOWER="$2"; shift 2;;
    -u) UPPER="$2"; shift 2;;
    -r) RADIAL=1; shift 1;;
    -d) DIPOLE=1; shift 1;;
    -e) EIGENF=1;SAVE=1; shift 1;;
    -f) FGONG=1; shift 1;;
    -s) SAVE=1; shift 1;;
    -R) RESOLUTION=1; shift 1;;
    -C) CONVERT=1; shift 1;;
    -U) UNITS="$2"; shift 2;;
    
    *) if [ -z "$INPUT" ]; then 
           INPUT="$1"
           shift 1
         else 
           echo "unknown option: $1" >&2
           exit 1
       fi;
  esac
done

if [ $HELP -gt 0 ] || [ -z "$INPUT" ]; then
    echo "Converter for .GYRE files to oscillation mode frequencies."
    echo "Usage: ./gyre2freqs.sh -i input -o output -t threads"
    echo "Flags: -s : save calculations directory"
    echo "       -e : calculate eigenfunctions (automatically turns on -s)"
    echo "       -r : only calculate radial modes"
    echo "       -d : only calculate dipole modes"
    echo "       -f : FGONG file format"
    echo "       -l : lower bound on frequency search"
    echo "       -u : upper bound on frequency search"
    echo "       -U : units such as 'CYC_PER_DAY' (default) or 'UHZ'"
    echo "       -R : increase grid resolution"
    exit
fi

if [[ $GYRE_DIR == '' ]]; then 
    echo "Error: GYRE_DIR not set!"
    exit 1
fi

## Check that the first input (GYRE file) exists
if [ ! -e "$INPUT" ]; then
    echo "Error: Cannot locate GYRE file $INPUT"
    exit 1
fi

## Pull out the name of the GYRE file
bname="$(basename $INPUT)"
fname="${bname%%.*}-freqs"
pname="${bname::-5}"

## If the OUTPUT argument doesn't exist, create a path from the filename 
if [ -z ${OUTPUT+x} ]; then
    path=$(dirname "$INPUT")/"$fname"
  else
    path="$OUTPUT"
fi

MODES="
&mode
    l=1
/
"
GRID_TYPE="INVERSE"

if [ $RESOLUTION -gt 0 ]; then 
    N_FREQ=1000 #10000
    GRID="
&grid
    w_ctr = 50
    w_osc = 50
    w_exp = 10
/
"
else
    N_FREQ=2000
    GRID="
&grid
    w_ctr = 10
    w_osc = 10
    w_exp = 2
/
"
fi

FORMAT="'MESA'"

MODE_ITEM_LIST=''
if [ $EIGENF -gt 0 ]; then
    SAVE=1
    MODE_ITEM_LIST="detail_file_format = 'TXT'
    detail_template = '%L_%N'
    detail_item_list = 'M_star,R_star,l,n_pg,n_p,n_g,freq,E,E_p,E_g,E_norm,M_r,x,xi_r,xi_h,dE_dx'"
fi

## Create a directory for the results and go there
mkdir -p "$path" 
cd "$path" 

logfile="gyre-l0.log"
#exec > $logfile 2>&1

## Create a gyre.in file to find the large frequency separation
echo "&model
    model_type = 'EVOL'
    file = '../$bname'
    file_format = $FORMAT
/

&constants
/

$MODES

&osc
    outer_bound = 'VACUUM'
    variables_set = 'JCD'
    inertia_norm = 'BOTH'
    x_ref = 1
/

&num
    diff_scheme = 'MAGNUS_GL6'
/

&scan
    grid_type = '$GRID_TYPE'
    freq_min_units = '$UNITS'
    freq_max_units = '$UNITS'
    freq_min = $LOWER
    freq_max = $UPPER
    n_freq = $N_FREQ
/

$GRID

&rot
/

&ad_output
    summary_file = '$fname.dat'
    summary_file_format = 'TXT'
    summary_item_list = 'l,m,n_pg,n_p,n_g,freq,E_norm'
    freq_units = '$UNITS'
    $MODE_ITEM_LIST
/

&nad_output
/

" >| "gyre.in"

# Run GYRE
$GYRE_DIR/bin/gyre gyre.in &>gyre.out

### Hooray!
if [ -f "$fname.dat" ]; then
    cp "$fname.dat" ..

    n_done=$(ls -l ../*-freqs.dat | wc -l)
    n_tot=$(ls -l ../*.data.GYRE | wc -l)
    echo "Conversion complete. Results in $fname.dat"

    if [ $SAVE -gt 0 ]; then exit 0; fi
    rm -rf *
    currdir=$(pwd)
    cd ..
    rm -rf "$currdir"
    exit 0
fi
echo "Something went wrong. Perhaps the following will help:"
cat gyre.out
cd ..
