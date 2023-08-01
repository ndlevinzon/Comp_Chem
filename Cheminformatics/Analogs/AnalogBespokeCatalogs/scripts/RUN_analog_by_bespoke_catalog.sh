#! /bin/bash

analog_rxns_json_file=$1 # Json file containing all the bespoke analoging reactions
analog_synthonRules_json_file=$2 # Json file containing all the synthon type inclusion/exclusion rules
input_mol_file=$3 # File containing list of all compounds to be analoged in form of: smi (\s or \t) cid
tc_max=$4 # Max Tc to find alternative BBs for

export SCRIPT_DIR='/nfs/soft/analog_by_bespoke_catalog/scripts'

# Split compounds in input mol file into their own files
while read line; do
	cid=$(echo $line | awk '{print $2}')
	echo $line > ${cid}_input
	echo ${cid}_input >> ${input_mol_file}_filenames # Save list of all input compound file names
	echo ${cid}_input >> made_tmp_files
	echo ${input_mol_file}_filenames >> made_tmp_files
done<$input_mol_file

# Source env with rdkit
source /mnt/nfs/ex9/work/ttummino/miniconda/etc/profile.d/conda.sh
conda activate base3.7

# Run Analoging Breakdown Rxns Script
for analog_mol_infile in $(cat ${input_mol_file}_filenames); do # Could imagine here is where you split into parallel jobs
	python3 $SCRIPT_DIR/Bespoke_Analoging_BreakdownRxns_1.0.py $analog_rxns_json_file $analog_mol_infile 
	# Output is a file named ${$analog_mol_infile}_analoging
	echo ${analog_mol_infile}_analoging >> breakdown_BB_filenames
	echo ${analog_mol_infile}_analoging >> made_tmp_files
	echo breakdown_BB_filenames >> made_tmp_files
done



# Source new env that can search Arthor (must be ssh'ed on gimel2)
source /nfs/home/ak87/exa/UCSF/SynthI/BESPOKE/arthor-env/bin/activate

for breakdown_BB_file in $(cat breakdown_BB_filenames); do
	while read line; do
		rxn_name=$(echo $line | awk '{print $2}')
		bbs_mol_list=$(echo $line | awk '{print $4}')
		bbs_synthonID_list=$(echo $line | awk '{print $5}')
		python3 $SCRIPT_DIR/analogs-arthor1_SVmod.py -i $bbs_mol_list -s $bbs_synthonID_list -t $tc_max -o ${breakdown_BB_file}_${rxn_name}_Tc${tc_max}_alternativeBBs -j $analog_synthonRules_json_file 

	done<$breakdown_BB_file
done


echo "Cleaning up tmp files..."
for tmp_file in $(cat made_tmp_files); do
	echo $tmp_file
	rm $tmp_file
done
rm made_tmp_files
