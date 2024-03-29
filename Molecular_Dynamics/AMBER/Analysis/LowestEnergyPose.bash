#!/bin/bash
#SBATCH --job-name=lowestE
#SBATCH --time=48:00:00
#SBATCH --mem=22000
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=owner-guest
#SBATCH --partition=lonepeak-guest
set -x

export AMBERHOME="/uufs/chpc.utah.edu/common/home/u0818159/amber/GIT/amber20"

perl mdout.perl ./mdout.1

# Extracts Lowest Energy Frame from summary.EPTOT
lowestEnergyFrame=$(cat summary.EPTOT | awk '{if($2<min) {min=$2;print $1"   "min}}' | tail -1 | awk '{print $1}')
echo "Lowest Energy Frame = $lowestEnergyFrame"

# Finds the Instance of Lowest Energy Frame in MDOUT file
line=$(grep $lowestEnergyFrame mdout.1)

# Parses NSTEP to Extract NSTEP Value and Divides by NTWX from MD Input File
if [[ ! -z "$line" ]]; then
    nstep=$(echo $line | awk '{print $3}')
    echo "NSTEP = $nstep"
else
    echo "Error: lowestEnergyFrame not found in mdout.1"
fi
nstep=$((nstep / 500))

cpptraj<<EOF
parm ../complex.topo
trajin ../TRAJ/traj.1 $nstep $nstep
trajout lowest_energy_struct.pdb pdb 
run
EOF

##############################################################
#______  _____   _          ______   ______ ______  ________ #
#|  ___ \(____ \ | |        (_____ \ / __   (_____ \(_______/#
#| |   | |_   \ \| |          ____) ) | //| | ____) )  ____  #
#| |   | | |   | | |         /_____/| |// | |/_____/  (___ \ #
#| |   | | |__/ /| |_____    _______|  /__| |_______ _____) )#
#|_|   |_|_____/ |_______)  (_______)\_____/(_______|______/ #
##############################################################                                                           
