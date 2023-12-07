#!/bin/bash

# Set the working directory to the current directory
MY_WORKDIR=`pwd`

# Temperature and Hamiltonian info.
TEMPERATURES=Temperatures.dat
HAMILTONIANS=Hamiltonians.dat
H_DIM=`cat $HAMILTONIANS | wc -l`
T_DIM=`cat $TEMPERATURES | wc -l`
echo "$T_DIM temperature dimensions."
echo "$H_DIM Hamiltonian dimensions."

# Check if dimensions are non-zero
if [[ $T_DIM -lt 1 || $H_DIM -lt 1 ]] ; then
  echo "Error: One or more dimensions are 0." > /dev/stderr
  exit 1
fi

# Read temperature and Hamiltonian variables
TVARS=`cat $TEMPERATURES`
HVARS=`cat $HAMILTONIANS`

# Starting coordinates location, expected file format is XXX.rst7
CRD_DIR=../starting_structures

# Total number of runs to execute
TOTAL_NRUN=10

# Create INPUT and groupfile for each run that will execute.
for (( RUN = 0; RUN < $TOTAL_NRUN; RUN++ )) ; do
  cd $MY_WORKDIR
  REXT=`printf "%03i" $RUN`
  RDIR="run.$REXT"
  echo "Creating run directory $RDIR"

  # Check if the run directory already exists
  if [[ -e $RDIR ]] ; then
    echo "Error: run directory exists." > /dev/stderr
    exit 1
  fi

  mkdir $RDIR
  cd $RDIR

  # Create groupfile and input for all replicas in this run
  if [[ -e groupfile ]] ; then 
    rm groupfile
  fi

  if [[ ! -e INPUT ]] ; then
    mkdir INPUT
  fi

  # Overall replica number
  REP=0   

  # Begin loop over Hamiltonians
  DIM=0

  for TOPFILE in $HVARS ; do
    TOP="../$TOPFILE"
    # Sanity check - make sure Top exists!
    if [[ ! -e $TOP ]] ; then
      echo "Error in Hamiltonian dim: $TOP does not exist." > /dev/stderr
      exit 1
    fi

    # Begin loop over Temperatures 
    for T in $TVARS ; do
      R=`printf "%03i" $REP`
      cat > INPUT/in.$R <<EOF
REMD for $TOPFILE at $T K ($REP)
&cntrl
    ntpr=1000, ntwr=750, ntwx=500, ntxo=2,
    ntf=2, ntc=2, ntp=0, ntb=1,
    ntt = 3, ig=-1, gamma_ln=2.0, cut = 8.0,
    nstlim=500, dt=0.002, numexchg=50000,
    imin=0, ntx=5, irest=1, ioutfm=1, iwrap=1,
    temp0=$T,   
/
EOF
      # Sanity check - if run 0 make sure input coords exist.
      if [[ $RUN -eq 0 && ! -e $CRD_DIR/$R.rst7 ]] ; then
        echo "Error: Input coords for run 0 '$CRD_DIR/$R.rst7' not found." > /dev/stderr
        echo "Error: Make sure the path $CRD_DIR is relative to subdirectory $RDIR" > /dev/stderr
        exit 1
      fi
      echo "-O -remlog rem.log -i INPUT/in.$R -p $TOP -c $CRD_DIR/$R.rst7 -o OUTPUT/rem.out.$R -inf INFO/reminfo.$R -r RST/$R.rst7 -x TRAJ/rem.crd.$R -l LOG/logfile.
$R" >> groupfile
      ((REP++))
    done # End T-loop
    ((DIM++))
  done # End H-loop

  # Since REP starts from one above, REP-1 is actual # REPS
  ((REP--))
  echo "$REP total replicas"

  # Generate remd.dim file
  cat > remd.dim <<EOF
Temperature REMD
&multirem
   exch_type = 'TEMPERATURE',
EOF
  #Create a temperature group for each Hamiltonian
  TGROUP=1
  for ((NH = 1; NH <= 8; NH++)) ; do
    printf "   group($NH,:) = " >> remd.dim
    for (( NT = 1; NT <= 24; NT++ )) ; do
      printf "$TGROUP," >> remd.dim
      ((TGROUP++))
    done
    printf "\n" >> remd.dim
  done
  cat >> remd.dim <<EOF
   desc = 'Temperature exchange from 300K to $T K'
/
Hamiltonian REMD
&multirem
   exch_type='HAMILTONIAN',
EOF
  # Create pairs of same temperature, different top
  for (( NT=1; NT <= 24; NT++ )) ; do
    printf "   group($NT,:) = " >> remd.dim
    for (( NH=$NT; NH <= $REP; NH += 24 )) ; do
      printf "$NH," >> remd.dim
    done
    printf "\n" >> remd.dim
  done
  cat >> remd.dim <<EOF
   desc = 'Scale vdWs'
/
EOF
  # Create output directories
  for DIR in OUTPUT TRAJ RST INFO LOG ; do
    mkdir $DIR
  done
  # Set striping for the TRAJ/OUTPUT directories
  echo "Setting striping for TRAJ and OUTPUT directories."
  lfs setstripe --count 160 TRAJ/
  lfs setstripe --count 160 OUTPUT/
  # Input coordinates for next run will be restarts of this
  CRD_DIR=../$RDIR/RST
done # End loop over runs
exit 0
