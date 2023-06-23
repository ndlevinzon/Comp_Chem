#!/bin/bash

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
		# handle booleans specially here- aim is to reduce user confusion
		if [ "$default" = "false" ] || [ "$default" = "true" ]; then # determine whether this is a boolean value based on the default
			export $env_name=$(to_lower ${!env_name}) # lowercase the boolean, just in case we have python enjoyers
			if ! [ ${!env_name} = "true" ] && ! [ ${!env_name} = "false" ]; then
				echo "expected boolean value for $env_name, must be true or false (not case sensitive). current value=${!env_name}"
				failed=1
			fi
		fi
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
exists_warning USE_PARALLEL_ARGS "addtl arguments for GNU parallel command" ""
exists_warning USE_CHARITY "use charity engine" "false"
exists_warning CHARITY_AUTHKEY "authentication key for charity engine" ""

#!~~QUEUE TEMPLATE~~!#
#exists_warning USE_MY_QUEUE "use MY_QUEUE" "false"
#exists_warning USE_MY_QUEUE_ARGS "addtl arguments for MY_QUEUE" ""

n=0
[ "$USE_SLURM" = "true" ] && n=$((n+1))
[ "$USE_SGE" = "true" ] && n=$((n+1))
[ "$USE_PARALLEL" = "true" ] && n=$((n+1))
[ "$USE_CHARITY" = "true" ] && n=$((n+1))

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
exists_warning SHRTCACHE "temporary local storage for job files" /dev/shm

echo "=================miscellaneous================="
exists_warning SUBMIT_WAIT_TIME "how many seconds to wait before submitting" 5
exists_warning USE_CACHED_SUBMIT_STATS "only check completion for jobs submitted in the latest iteration. Faster re-submission, but will ignore jobs that have been manually reset" "false"

if [ $failed -eq 1 ]; then
	echo "exiting with error"
	exit 1
fi

mkdir -p $EXPORT_DEST/logs
n=1
njobs=0

function handle_input_source {
	exts=$@
	if [ -d $INPUT_SOURCE ]; then
		printf "" > $EXPORT_DEST/file_list
		for ext in $exts; do
			find $INPUT_SOURCE -name '*.'"$ext" -type f | sort >> $EXPORT_DEST/file_list
		done
	else
		cp $INPUT_SOURCE $EXPORT_DEST/file_list
	fi
	echo $EXPORT_DEST/file_list
}

if [ "$USE_DB2_TGZ" = "true" ]; then
	inp=$(handle_input_source 'db2.tgz' 'db2.tar.gz')
	n_input_tot=$(cat $inp | wc -l)
	get_input_cmd="seq 1 $(( (n_input_tot+USE_DB2_TGZ_BATCH_SIZE-1)/USE_DB2_TGZ_BATCH_SIZE))"
elif [ "$USE_DB2" = "true" ]; then
	inp=$(handle_input_source 'db2.gz' 'db2')
	n_input_tot=$(cat $inp | wc -l)
	get_input_cmd="seq 1 $(( (n_input_tot+USE_DB2_BATCH_SIZE-1)/USE_DB2_BATCH_SIZE))"
else
	echo "you must select an input type! set USE_DB2=true or USE_DB2_TGZ=true"
	exit 1
fi

# new variable- keep track of each resubmission for posterity, recorded in joblist files
# find this value by counting joblist files
RESUBMIT_COUNT=0
while [ -f $EXPORT_DEST/joblist.$RESUBMIT_COUNT ]; do
	RESUBMIT_COUNT=$((RESUBMIT_COUNT+1))
done
# resubmitting for the 2nd time onward will be faster if we use the cached joblist from the previous run. optional feature 
if [ $RESUBMIT_COUNT -gt 0 ] && [ $USE_CACHED_SUBMIT_STATS = "true" ]; then
	get_input_cmd="cat $EXPORT_DEST/joblist.$((RESUBMIT_COUNT-1))"
fi

nrestart=0
for input in $($get_input_cmd); do
	if [ $RESUBMIT_COUNT -eq 0 ]; then # first run don't bother checking file existence
		echo $input
		njobs=$((njobs+1))
	elif ! [ -f $EXPORT_DEST/$input/OUTDOCK.0 ]; then
		echo $input
		njobs=$((njobs+1))
	elif ! [ -f $EXPORT_DEST/$input/test.mol2.gz.0 ]; then
		#rm $EXPORT_DEST/$input/*
		echo $input
		njobs=$((njobs+1))
	elif [ -f $EXPORT_DEST/$input/OUTDOCK.0 ] && [ -f $EXPORT_DEST/$input/restart ]; then
		echo $input
		njobs=$((njobs+1))
		nrestart=$((nrestart+1))
	fi
done > $EXPORT_DEST/joblist.$RESUBMIT_COUNT

echo "attempt number $RESUBMIT_COUNT:"
echo "submitting $njobs out of $input jobs over $n_input_tot files. $((input-njobs)) already complete, $((nrestart)) partially complete"

[ -z $QSUB_EXEC ] && QSUB_EXEC=qsub
[ -z $SBATCH_EXEC ] && SBATCH_EXEC=sbatch
[ -z $PARALLEL_EXEC ] && PARALLEL_EXEC=parallel
[ -z $CHARITY_EXEC ] && CHARITY_EXEC=ce-cli

#!~~QUEUE TEMPLATE~~!#
#[ -z $MY_QUEUE_EXEC ] && MY_QUEUE_EXEC=something

echo "!!! save the following to its own file for a re-usable superscript !!!"
echo "==============================================================="
echo "#!/bin/bash"
for var in EXPORT_DEST INPUT_SOURCE DOCKFILES DOCKEXEC \
 SHRTCACHE SHRTCACHE_USE_ENV \
 USE_DB2_TGZ USE_DB2_TGZ_BATCH_SIZE USE_DB2 USE_DB2_BATCH_SIZE\
 USE_SLURM USE_SLURM_ARGS USE_SGE USE_SGE_ARGS USE_PARALLEL USE_PARALLEL_ARGS MAX_PARALLEL\
 QSUB_EXEC SBATCH_EXEC PARALLEL_EXEC \
 SUBMIT_WAIT_TIME USE_CACHED_SUBMIT_STATS; do
	echo "export $var=${!var}"
done
echo "bash $BINPATH"
echo "==============================================================="

SGE_ENV_ARGS=""
SLURM_ENV_ARGS=""
# pass in rundock specific vars here- passing in the subdock specific ones as well runs into issue, mostly because of USE_SGE_ARGS and the like having non-standard formatting (damn whitespace...)
for var in EXPORT_DEST INPUT_SOURCE DOCKFILES DOCKEXEC \
 SHRTCACHE SHRTCACHE_USE_ENV \
 USE_DB2_TGZ USE_DB2_TGZ_BATCH_SIZE USE_DB2 USE_DB2_BATCH_SIZE \
 USE_SLURM USE_SGE USE_PARALLEL \
 RESUBMIT_COUNT; do
	#!~~QUEUE TEMPLATE~~!#
        # your queueing system may require explicit enumeration of environment values to export (like sge)
        # add a similar implementation here if required
 	export $var
	[ -z "$SGE_ENV_ARGS" ] && SGE_ENV_ARGS="-v $var=${!var}" || SGE_ENV_ARGS="$SGE_ENV_ARGS -v $var=${!var}"
	# explicitly define slurm args here- don't copy whole user environment with --export=ALL, in case something weird is there
	[ -z "$SLURM_ENV_ARGS" ] && SLURM_ENV_ARGS="--export=$var" || SLURM_ENV_ARGS="$SLURM_ENV_ARGS,$var"
done


if [ $njobs -eq 0 ]; then
	echo "all $input jobs complete!"
	rm $EXPORT_DEST/joblist.$RESUBMIT_COUNT
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
USE_PARALLEL_ARGS="$USE_PARALLEL_ARGS --termseq USR1,10000,KILL,25" # USR1, wait 10s, then kill, then 25ms and exit
#USE_PARALLEL_ARGS="$USE_PARALLEL_ARGS --termseq TERM,10000,KILL,25" # uncomment for testing unexpected interrupts

if [ "$USE_SGE" = "true" ]; then
	[ $MAX_PARALLEL -gt 0 ] && USE_SGE_ARGS="$USE_SGE_ARGS -tc $MAX_PARALLEL"
echo	$QSUB_EXEC $SGE_ENV_ARGS $SGE_LOG_ARGS -cwd -S /bin/bash -t 1-$njobs $USE_SGE_ARGS $RUNDOCK_PATH
	$QSUB_EXEC $SGE_ENV_ARGS $SGE_LOG_ARGS -cwd -S /bin/bash -t 1-$njobs $USE_SGE_ARGS $RUNDOCK_PATH

elif [ "$USE_PARALLEL" = "true" ]; then
	[ $MAX_PARALLEL -gt 0 ] && USE_PARALLEL_ARGS="$USE_PARALLEL_ARGS -j $MAX_PARALLEL"
	export JOB_ID='test'
echo	env time -v -o $EXPORT_DEST/perfstats $PARALLEL_EXEC --results $EXPORT_DEST/logs $USE_PARALLEL_ARGS bash $RUNDOCK_PATH ::: $(seq 1 $njobs)
	env time -v -o $EXPORT_DEST/perfstats $PARALLEL_EXEC --results $EXPORT_DEST/logs $USE_PARALLEL_ARGS bash $RUNDOCK_PATH ::: $(seq 1 $njobs)

elif [ "$USE_SLURM" = "true" ]; then
	USE_SLURM_ARGS="$USE_SLURM_ARGS --array=1-$njobs"
	[ $MAX_PARALLEL -gt 0 ] && USE_SLURM_ARGS="${USE_SLURM_ARGS}%$MAX_PARALLEL"
		#                            vv causes slurm to interrupt script 10s prior to timeout (if one is specified)
echo	$SBATCH_EXEC $USE_SLURM_ARGS --signal=B:USR1@10 $SLURM_ENV_ARGS $USE_SLURM_ARGS $SLURM_LOG_ARGS $RUNDOCK_PATH
	$SBATCH_EXEC $USE_SLURM_ARGS --signal=B:USR1@10 $SLURM_ENV_ARGS $USE_SLURM_ARGS $SLURM_LOG_ARGS $RUNDOCK_PATH

elif [ "$USE_CHARITY" = "true" ]; then
	[ $MAX_PARALLEL -lt 1  ] && MAX_PARALLEL=1  && echo MAX_PARALLEL too small, defaulting to MAX_PARALLEL=1
	[ $MAX_PARALLEL -gt 50 ] && MAX_PARALLEL=50 && echo MAX_PARALLEL too large, defaulting to MAX_PARALLEL=50
	# check for errors etc.
	if [ -z $CHARITY_AUTHKEY ]; then
		echo "need an authentication key for charity!"
		exit 1
	elif ! [ $USE_DB2_TGZ = "true" ]; then
		echo "must select USE_DB2_TGZ!"
		exit 1
	elif ! [ $USE_DB2_TGZ_BATCH_SIZE = 1 ]; then
		echo "USE_DB2_TGZ_BATCH_SIZE must be 1!"
		exit 1
	fi
	if ! [[ $DOCKFILES == *.tgz ]]; then
		if ! [ -f ${DOCKFILES}.tgz ]; then
			echo "creating ${DOCKFILES}.tgz"
			pushd $(dirname $DOCKFILES) 1>&2 2>/dev/null
			tar -czf ${DOCKFILES}.tgz $(basename $DOCKFILES)
			popd 1>&2 2>/dev/null
		fi
                DOCKFILES=${DOCKFILES}.tgz
        fi
echo    $PARALLEL_EXEC -j $MAX_PARALLEL -a $EXPORT_DEST/joblist.$RESUBMIT_COUNT -a $EXPORT_DEST/file_list $CHARITY_EXEC --debug true --app docker:dockingorg/dock_ce:latest --env RESUBMIT_COUNT=$RESUBMIT_COUNT USE_CHARITY=true --commandline "bash /bin/rundock.bash" --auth $CHARITY_AUTHKEY --inputfile $DOCKFILES {2} --outputdir $EXPORT_DEST/{1}
        $PARALLEL_EXEC -j $MAX_PARALLEL -a $EXPORT_DEST/joblist.$RESUBMIT_COUNT -a $EXPORT_DEST/file_list $CHARITY_EXEC --debug true --app docker:dockingorg/dock_ce:latest --env RESUBMIT_COUNT=$RESUBMIT_COUNT USE_CHARITY=true --commandline "bash /bin/rundock.bash" --auth $CHARITY_AUTHKEY --inputfile $DOCKFILES {2} --outputdir $EXPORT_DEST/{1}
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
	
