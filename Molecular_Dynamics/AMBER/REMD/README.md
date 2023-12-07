If you’re reading this, you are trying to do a multidimensional replica exchange (likely on Frontera). 

First, you will need to determine your model system and which Hamiltonian (dihedral, vdW, etc) you will change. 

Build your systems (topologies and coordinates) as you usually use tleap. Be sure to note the number of solute atoms and solution molecules. You will likely use parmed to manipulate your original topology into the other topologies associated with your other Hamiltonians. You will need an even number of Hamiltonians. Create a “Hamiltonian.dat” file and list every topology, each on a separate line. For example:
10.topo 
7.topo 
etc…

Then, make your way to this website: https://virtualchemistry.org/remd-temperature-generator/. Here you will calculate the temperatures. Your exchange probability will be 0.2. Choose an upper and lower limit (for RNA I did 270-300). Add your water number and atoms in your protein (or nucleic acid). Then, hit calculate. You will be given a list of temperatures. This should be an even number. Take these temperatures and create a “Temperature.dat” file where each temperature is listed, each on a separate line. For example:
270.3
272.9
Etc…

Make a directory called “starting_structures” and minimize/equilibrate the system as you would normally for each Hamiltonian/temperature pair. For example, if you have 8 Hamiltonians and 24 temperatures you will get 192 combinations. Make sure the final restart files are labeled starting at 000.rst7 and going all the way up to the one less than the total number of combinations (bc you started at 0). For this example, mine would end at 191.rst7 .

Next, using the bash script provided below we are going to make all of the directories and groupfiles and input files necessary to run. 


#!/bin/bash
MY_WORKDIR=`pwd`
# Temperature and Hamiltonian info.
TEMPERATURES=Temperatures.dat
HAMILTONIANS=Hamiltonians.dat
H_DIM=`cat $HAMILTONIANS | wc -l`
T_DIM=`cat $TEMPERATURES | wc -l`
echo "$T_DIM temperature dimensions."
echo "$H_DIM Hamiltonian dimensions."
if [[ $T_DIM -lt 1 || $H_DIM -lt 1 ]] ; then
  echo "Error: One or more dimensions are 0." > /dev/stderr
  exit 1
fi
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

This should create 10 running directories. Within the running directory you will have INPUT/ OUTPUT/ RST/ TRAJ/ INFO/ and LOG/ directories. You will also have a groupfile and a remd.dim. The remd.dim file is how your temperatures and hamiltonians get grouped. I had issues with my last system (192) being left off of the 24th group of the Hamiltonian grouping. This will cause an error so you may want to double-check that all the groups have the same number of systems. 

Now all you need is to run the simulations! 

This will require “pmemd.cuda.MPI” if you are running on GPUs or “pmemd.MPI” for running on CPUs. If you have the computer availability, you should run on GPUs bc it will be sooooo much faster. I have Amber22 compiled on Frontera and CHPC has Amber 20 as a module you can use. I will give example scripts for the following

 
For running on Frontera you will want to submit this script:

#!/bin/bash
#SBATCH -J Mremd-OL3
#SBATCH -p rtx
#SBATCH -N 16
#SBATCH -n 192
#SBATCH -t 24:00:00
#SBATCH -A MCB20008
#SBATCH --mail-type=all
#SBATCH --mail-user=lauren.winkler@utah.edu
set -x
module load mvapich2-x/2.3
module load cuda/11.0

DIR="OUT"
ibrun -np 192 /home1/07450/lgw_19/amber22/amber22/bin/pmemd.cuda.MPI -ng 192  -groupfile groupfile -remd-file remd.dim -remlog remlog

OPTIONAL, but helpful:
#SBATCH --dependency=5639832

Frontera does not let you submit the next script from your current script (no “sbatch prod.bash” at the end of the script like we can do at CHPC). So I would submit them all at once using the dependency feature so that when one stopped, the next one was already queued up. 


Also you will want to run this in one of your “scratch” directories. Remember that scratch directories are wiped about once a week if they aren’t touched, according to TACC. (I find this isn’t always true, maybe once every couple of months) Regardless, you should be diligent about “touching” your files regularly or transferring your directories over to CHPC for safekeeping. I would recommend checking out the data transfer nodes for easy transferring (https://chpc.utah.edu/documentation/software/slurm-datatransfernode.php#:~:text=CHPC%20now%20has%20enabled%20the%20use%20of%20the,RAM%20are%20made%20available%20to%20run%20SLURM%20jobs.). 

For running on CHPC- You likely will not be able to run on CHPC bc of the limited number of available GPUs but I will provide a sample script in case you want to do HREMD or TREMD, which require less computer power but takes longer to run. 

#!/bin/bash
#SBATCH --job-name=hremd-dihedral
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --tasks=8
#SBATCH --mem=0
#SBATCH --gres=gpu:2080ti:1
#SBATCH --account=cheatham-gpu-np
#SBATCH --partition=cheatham-gpu-np
set -x
module load gcc/8.5.0
module load mvapich2/2.3.7
module load intel-oneapi-mpi/2021.4.0
module load cuda/11.0
module load amber/20.20-gpu
mpirun -np 8 pmemd.cuda.MPI -rem 3 -ng 8 -groupfile groupfile
cd ../run.102/
sbatch prod.bash


Another optional but helpful script- this one will write all the production scripts into all of your running directories, helpful when you have many directories like in T-REMD/H-REMD. This is a c-shell script. 

#!/bin/csh
set cnt = 0
set cntmax = 9
while ( ${cnt} <= ${cntmax} )
    @ ncnt = ${cnt} + 1
    cd run.00${cnt}
    cat > prod.bash <<EOF
#!/bin/bash
#SBATCH -J Mremd-OL3
#SBATCH -p rtx
#SBATCH -N 16
#SBATCH -n 192
#SBATCH -t 24:00:00
#SBATCH -A MCB20008
#SBATCH --mail-type=all
#SBATCH --mail-user=lauren.winkler@utah.edu
set -x
module load mvapich2-x/2.3
module load cuda/11.0

DIR="OUT"
ibrun -np 192 /home1/07450/lgw_19/amber22/amber22/bin/pmemd.cuda.MPI -ng 192  -groupfile groupfile -remd-file remd.dim -remlog remlog
EOF
    cd ../
@ cnt += 1
end


Once you’ve collected enough data, you will want to show that you have reached convergence. A quick way to do this is to show that each trajectory (unsorted) is accessing the same structures, shown by overlapping RMSD histograms for each trajectory. 

So for each run I was output 192 trajectories. I want to take all of the trajectories of the same number and make sure they sample the same space. In this example I had 8 runs,

#!/bin/bash
#SBATCH --job-name=remd-water
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --tasks=4
#SBATCH --mem=0
#SBATCH --gres=gpu:2080ti:1
#SBATCH --account=notchpeak-gpu
#SBATCH --partition=notchpeak-gpu
set -x
cpptraj='/uufs/chpc.utah.edu/common/home/u0818159/amber18/bin/cpptraj'
for n in {0..9}
do
cpptraj<<EOF 
parm ../dna.topo
trajin ../run.000/TRAJ/rem.crd.00$n
trajin ../run.001/TRAJ/rem.crd.00$n
trajin ../run.002/TRAJ/rem.crd.00$n
trajin ../run.003/TRAJ/rem.crd.00$n
trajin ../run.004/TRAJ/rem.crd.00$n
trajin ../run.005/TRAJ/rem.crd.00$n
trajin ../run.006/TRAJ/rem.crd.00$n
trajin ../run.007/TRAJ/rem.crd.00$n
trajin ../run.008/TRAJ/rem.crd.00$n
box nobox
trajout remd.reptraj.00$n nobox
go 
quit 
EOF
    
$AMBERHOME/bin/cpptraj <<EOF
parm ../dna.topo
trajin remd.reptraj.00$n 1 last 10
strip :Na+,WAT,Cl-
parm /uufs/chpc.utah.edu/common/home/cheatham-group4/lauren/metastable-DRUDE/GACC/GACC-rep-structures/AformMaj.pdb [Aform]
reference /uufs/chpc.utah.edu/common/home/cheatham-group4/lauren/metastable-DRUDE/GACC/GACC-rep-structures/AformMaj.pdb parm [Aform]
rms RNA :1-4&!@H= reference :1-4&!@H= mass out RMSD/rmsd00$n.dat
hist RNA norm min 0 max 20 bins 200 out RMSD/rmsd00$n.hist
go 
quit
EOF
Done

Plotting the resulting histograms on the same graph should show you all the replicas are sampling the same space (only showing 10 replicas here but you get the point)



Once you’ve established convergence within a run, you will want to establish convergence between two independent runs. 

Now you may want to look at specific properties. This next script will sort all the trajectories by hamiltonian and temperature pair. For example, the resulting remd.crd.0 will be the [1,1] pair. Aka the lowest hamiltonian and lowest temperature. In this example it will be all the time the simulation spend at the Hamiltonian “10.topo” and temperature 270K. Then remd.crd.1 would be [1,2]. This is the time the simulation spend at the Hamiltonian “7.topo” and temperature 270K, etc. 

#!/bin/bash
#SBATCH --job-name=remd-water
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --tasks=8
#SBATCH --mem=0
#SBATCH --gres=gpu:1080ti:1
#SBATCH --account=lonepeak-gpu
#SBATCH --partition=lonepeak-gpu
set -x
cpptraj='/uufs/chpc.utah.edu/common/home/u0950975/cpptraj-6.18.0'
cpptraj<<EOF 
parm ../dna.topo
ensemble ../run.000/TRAJ/rem.crd.000
ensemble ../run.001/TRAJ/rem.crd.000 
ensemble ../run.002/TRAJ/rem.crd.000 
ensemble ../run.003/TRAJ/rem.crd.000 
ensemble ../run.004/TRAJ/rem.crd.000 
ensemble ../run.005/TRAJ/rem.crd.000 
ensemble ../run.006/TRAJ/rem.crd.000 
ensemble ../run.007/TRAJ/rem.crd.000 
ensemble ../run.008/TRAJ/rem.crd.000 
strip :WAT
trajout by-H/remd.crd
go 
quit 
EOF

Clustering-
cpptraj<<EOF
parm ../../../1/dnanowat.topo
trajin ../../../1/analysis/by-H/remd.crd.16 1 last 10
trajin ../../../2/analysis/by-H/remd.crd.16 1 last 10
strip :WAT
cluster dbscan minpoints 25 epsilon 0.5 sievetoframe \
 :1@N2,O6,C1',P,:2@H2,N6,C1',P,:3@O2,H5,C1',P,:4@O2,H5,C1',P \
 sievetoframe sieve 3 out cvt.dat summary summary.dat \
    info info.dat \
  cpopvtime cpopvtime.agr normframe \
  repout rep repfmt pdb \
  singlerepout singlerep.nc singlerepfmt netcdf \
clusterout traj clusterfmt netcrd \
avgout avg avgfmt pdb
run
quit 

Good luck and let me know if you need any help or clarification!
