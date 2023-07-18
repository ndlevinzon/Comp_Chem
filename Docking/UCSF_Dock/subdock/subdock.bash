#!/bin/bash
# req:
# EXPORT_DEST
# INPUT_SOURCE
# DOCKFILES
# DOCKEXEC
# SHRTCACHE
# LONGCACHE

BINDIR=$(dirname $0)
BINDIR=${BINDIR-.}
BINPATH=$0
if ! [[ "$BINDIR" == /* ]]; then
	BINDIR=$PWD/$BINDIR
	BINPATH=$BINDIR/$(basename $0)
fi
RUNDOCK_PATH=$BINDIR/rundock.bash

function to_upper {
        echo "$1" | tr '[:lower:]' '[:upper:]' | tr '-' '_'
}
function to_lower {
        echo "$1" | tr '[:upper:]' '[:lower:]' | tr '_' '-'
}

failed=0
function exists {
	env_name=$1
	desc=$2
	if [ -z "${!env_name}" ]; then
		echo "expected env arg: $env_name, --$(to_lower $env_name)"
		echo "arg description: $desc"
		echo 
		failed=1
	else
		echo $env_name=${!env_name}
		echo
	fi
}

function exists_warning {
	env_name=$1
	desc=$2
	default=$3
	if [ -z "${!env_name}" ]; then
		echo "optional env arg missing: $env_name, --$(to_lower $env_name)"
		echo "arg description: $desc"
		echo "defaulting to $default"
		export $env_name="$default"
		echo
	else
		echo $env_name=${!env_name}
		echo
	fi
}

# allow arguments to be passed in this way, cause why not. live your life
# quotation marks around the $@ very important!!! here is why:
# cmd.bash a="b c" d
# $1 will be "a=b c", $2 will be "d"
# however, if you iterate through the arguments with $@, something odd happens
# $1 becomes "a=b", $2 becomes "c" and $3 becomes "d"
# unless you include the quotation marks around $@, in which case it reverts back to the "normal" behavior
for arg in "$@"; do
	# remove leading "--"
	arg_s=$(echo "$arg" | tail -c+3)
	var=$(echo $arg_s | cut -d'=' -f1)
	val=$(echo $arg_s | cut -d'=' -f2-)
	var=$(to_upper $var)
	export $var="$val"
done

echo "SUBDOCK! Run docking workloads via job controller of your choice"

echo "=================required arguments================="
exists EXPORT_DEST "nfs output destination for OUTDOCK and test.mol2.gz files"
exists INPUT_SOURCE "nfs directory containing one or more .db2.tgz files OR a file containing a list of db2.tgz files"
exists DOCKFILES "nfs directory containing dock related files and INDOCK configuration for docking run"
exists DOCKEXEC "nfs path to dock executable"

echo "=================job controller settings================="
# queue system is active if set to "true" otherwise inactive
exists_warning USE_SLURM "use slurm" "false"
exists_warning USE_SLURM_ARGS "addtl arguments for SLURM sbatch command" ""
exists_warning USE_SGE "use sge" "false"
exists_warning USE_SGE_ARGS "addtl arguments for SGE qsub command" ""
exists_warning USE_PARALLEL "use GNU parallel" "false"

#!~~QUEUE TEMPLATE~~!#
#exists_warning USE_MY_QUEUE "use MY_QUEUE" "false"
#exists_warning USE_MY_QUEUE_ARGS "addtl arguments for MY_QUEUE" ""

n=0
[ "$USE_SLURM" = "true" ] && n=$((n+1))
[ "$USE_SGE" = "true" ] && n=$((n+1))
[ "$USE_PARALLEL" = "true" ] && n=$((n+1))

#!~~QUEUE TEMPLATE~~!#
#[ "$USE_MY_QUEUE" = "true" ] && n=$((n+1))

if [ $n -gt 1 ]; then
	echo "cannnot select more than one job controller!"
	exit 1
elif [ $n -lt 1 ]; then
	echo "must select a job controller!"
	exit 1
fi

echo "=================input settings================="
exists_warning USE_DB2_TGZ "dock db2.tgz tar files" "true"
exists_warning USE_DB2_TGZ_BATCH_SIZE "how many db2.tgz to evaluate per batch" 1
exists_warning USE_DB2 "dock db2.gz individual files" "false"
exists_warning USE_DB2_BATCH_SIZE "how many db2.gz to evaluate per batch" 100

if [ "$USE_DB2_TGZ" = "true" ] && [ "$USE_DB2" = "true" ]; then
	echo "cannot select more than one input type!"
	exit 1
elif ! [ "$USE_DB2_TGZ" = "true" ] && ! [ "$USE_DB2" = "true" ]; then
	echo "must select an input type!"
	exit 1
fi

echo "=================addtl job configuration================="
exists_warning MAX_PARALLEL "max jobs allowed to run in parallel" "-1"
exists_warning SHRTCACHE "temporary local storage for job files" /scratch
exists_warning LONGCACHE "longer term storage for files shared between jobs" /scratch

echo "=================miscellaneous================="
exists_warning SUBMIT_WAIT_TIME "how many seconds to wait before submitting" 5

if [ $failed -eq 1 ]; then
	echo "exiting with error"
	exit 1
fi

# old dummy-proofing code, not needed anymore as we programmatically fix this mistake now
#ligand_atom_file=$(grep -w ligand_atom_file $DOCKFILES/INDOCK | awk '{print $2}')
#if [ "$ligand_atom_file" != "split_database_index" ]; then
#	echo "ligand_atom_file option in INDOCK is [$ligand_atom_file], should be [split_database_index]. Please change and try again."
#	exit 1
#fi

if [ -w $DOCKFILES ]; then
	cat $DOCKFILES/* | sha1sum | awk '{print $1}' > $DOCKFILES/.shasum
fi

mkdir -p $EXPORT_DEST/logs
n=1
njobs=0
warned=

function handle_input_source {
	exts=$1
	if [ -d $INPUT_SOURCE ]; then
		printf "" > $EXPORT_DEST/file_list
		for ext in exts; do
			find $INPUT_SOURCE -name "$ext" -type f | sort >> $EXPORT_DEST/file_list
		done
	else
		cp $INPUT_SOURCE $EXPORT_DEST/file_list
	fi
	echo $EXPORT_DEST/file_list
}

if [ "$USE_DB2_TGZ" = "true" ]; then
	inp=$(handle_input_source '*.db2.tgz')
	n_input_tot=$(cat $inp | wc -l)
	get_input_cmd="seq 1 $(( (n_input_tot+USE_DB2_TGZ_BATCH_SIZE-1)/USE_DB2_TGZ_BATCH_SIZE))"
elif [ "$USE_DB2" = "true" ]; then
	inp=$(handle_input_source '*.db2.gz')
	n_input_tot=$(cat $inp | wc -l)
	get_input_cmd="seq 1 $(( (n_input_tot+USE_DB2_BATCH_SIZE-1)/USE_DB2_BATCH_SIZE))"
else
	echo "you must select an input type! set USE_DB2=true or USE_DB2_TGZ=true"
	exit 1
fi

for input in $($get_input_cmd); do
	if ! [ -f $EXPORT_DEST/$n/OUTDOCK.0 ]; then
		echo $input $n
		njobs=$((njobs+1))
	elif ! [ -f $EXPORT_DEST/$n/test.mol2.gz.0 ]; then
		rm $EXPORT_DEST/$n/*
		echo $input $n
		njobs=$((njobs+1))
	elif [ -f $EXPORT_DEST/$n/OUTDOCK.0 ] && [ -f $EXPORT_DEST/$n/restart ]; then
		echo $input $n
		njobs=$((njobs+1))
	fi
	n=$((n+1))
done > $EXPORT_DEST/joblist
n=$((n-1))

echo "submitting $njobs out of $n jobs over $n_input_tot files. $((n-njobs)) already complete"

[ -z $QSUB_EXEC ] && QSUB_EXEC=qsub
[ -z $SBATCH_EXEC ] && SBATCH_EXEC=sbatch
[ -z $PARALLEL_EXEC ] && PARALLEL_EXEC=parallel

#!~~QUEUE TEMPLATE~~!#
#[ -z $MY_QUEUE_EXEC ] && MY_QUEUE_EXEC=something


var_args=
echo "!!! save the following to its own file for a re-usable superscript !!!"
echo "==============================================================="
echo "#!/bin/bash"
for var in EXPORT_DEST INPUT_SOURCE DOCKFILES\
 DOCKEXEC SHRTCACHE LONGCACHE SHRTCACHE_USE_ENV\
 USE_DB2_TGZ USE_DB2_TGZ_BATCH_SIZE USE_DB2 USE_DB2_BATCH_SIZE\
 USE_SLURM USE_SLURM_ARGS USE_SGE USE_SGE_ARGS USE_PARALLEL MAX_PARALLEL\
 QSUB_EXEC SBATCH_EXEC PARALLEL_EXEC; do
	export $var
	echo "export $var=${!var}"
	
	#!~~QUEUE TEMPLATE~~!#
	# your queueing system may require explicit enumeration of environment values to export (like sge)
	# add a similar implementation here if required
done
# pass in rundock specific vars here- passing in the subdock specific ones as well runs into issue, mostly because of USE_SGE_ARGS and the like having non-standard formatting
for var in EXPORT_DEST INPUT_SOURCE DOCKFILES DOCKEXEC \
 SHRTCACHE LONGCACHE SHRTCACHE_USE_ENV \
 USE_DB2_TGZ USE_DB2_TGZ_BATCH_SIZE USE_DB2 USE_DB2_BATCH_SIZE \
 USE_SLURM USE_SGE USE_PARALLEL; do
	[ -z "$var_args" ] && var_args="-v $var=${!var}" || var_args="$var_args -v $var=${!var}"
done
echo "bash $BINPATH"
echo "==============================================================="


if [ $njobs -eq 0 ]; then
	echo "all $n jobs complete!"
	exit 0
fi

echo "sleeping for $SUBMIT_WAIT_TIME seconds before submitting, ctrl-C to stop"
if ! [ -z $SLURM_SETTINGS ]; then
        source $SLURM_SETTINGS
fi

sleep $SUBMIT_WAIT_TIME 
echo "submitting jobs"

SLURM_LOG_ARGS="-o $EXPORT_DEST/logs/%a.out -e $EXPORT_DEST/logs/%a.err"
SGE_LOG_ARGS="-o $EXPORT_DEST/logs -e $EXPORT_DEST/logs"

if [ "$USE_SGE" = "true" ]; then
	if [ $MAX_PARALLEL -gt 0 ]; then
	echo	$QSUB_EXEC $var_args $SGE_LOG_ARGS -cwd -S /bin/bash -t 1-$njobs -tc $MAX_PARALLEL $USE_SGE_ARGS $RUNDOCK_PATH
		$QSUB_EXEC $var_args $SGE_LOG_ARGS -cwd -S /bin/bash -t 1-$njobs -tc $MAX_PARALLEL $USE_SGE_ARGS $RUNDOCK_PATH
	else
	echo	$QSUB_EXEC $var_args $SGE_LOG_ARGS -cwd -S /bin/bash -t 1-$njobs $USE_SGE_ARGS $RUNDOCK_PATH
		$QSUB_EXEC $var_args $SGE_LOG_ARGS -cwd -S /bin/bash -t 1-$njobs $USE_SGE_ARGS $RUNDOCK_PATH
	fi
elif [ "$USE_PARALLEL" = "true" ]; then
	export JOB_ID='test'
	if [ $MAX_PARALLEL -gt 0 ]; then
		/usr/bin/time -v -o $EXPORT_DEST/perfstats $PARALLEL_EXEC --results $EXPORT_DEST/logs -j $MAX_PARALLEL bash $RUNDOCK_PATH ::: $(seq 1 $njobs)
	else
		/usr/bin/time -v -o $EXPORT_DEST/perfstats $PARALLEL_EXEC --results $EXPORT_DEST/logs bash $RUNDOCK_PATH ::: $(seq 1 $njobs)
	fi
elif [ "$USE_SLURM" = "true" ]; then
	if [ $MAX_PARALLEL -gt 0 ]; then
		#                            VV causes slurm to interrupt script 2m prior to timeout (if one is specified)
		$SBATCH_EXEC $USE_SLURM_ARGS --signal=B:USR1@120 --array=1-$njobs%$MAX_PARALLEL $SLURM_LOG_ARGS $RUNDOCK_PATH
	else
		$SBATCH_EXEC $USE_SLURM_ARGS --signal=B:USR1@120 --array=1-$njobs $SLURM_LOG_ARGS $RUNDOCK_PATH
	fi
#!~~QUEUE TEMPLATE~~!#
#elif [ "$USE_MY_QUEUE" = "true" ]; then
#	if [ $MAX_PARALLEL -gt 0 ]; then 
# 		# ntasks, script, etc. are dummy arguments- replace appropriately
#		$MY_QUEUE_EXEC $USE_MY_QUEUE_ARGS --ntasks=$njobs --max-parallel-tasks=$MAX_PARALLEL --script=$RUNDOCK_PATH
#	else
#		$MY_QUEUE_EXEC $USE_MY_QUEUE_ARGS --ntasks=$njobs --script=$RUNDOCK_PATH
#	fi
else
	echo "you need to specify at least one job submission method!!"
	exit 1
fi

exit 0
	
