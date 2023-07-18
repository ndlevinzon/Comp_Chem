#!/bin/bash
# req:
# EXPORT_DEST
# INPUT_SOURCE
# DOCKFILES
# DOCKEXEC
# SHRTCACHE
# LONGCACHE
# JOB_ID
# TASK_ID


function log {
	echo "[$(date +%X)]: $@"
}

if [ -z $SHRTCACHE_USE_ENV ]; then
	SHRTCACHE=${SHRTCACHE-/dev/shm}
else
	SHRTCACHE=${!SHRTCACHE_USE_ENV}
fi
LONGCACHE=${LONGCACHE-/scratch}

if [ "$USE_PARALLEL" = "true" ]; then
	TASK_ID=$1
elif [ "$USE_SLURM" = "true" ]; then
	JOB_ID=${SLURM_ARRAY_JOB_ID}
	TASK_ID=${SLURM_ARRAY_TASK_ID}
elif [ "$USE_SGE" = "true" ]; then
	JOB_ID=$JOB_ID
	TASK_ID=$SGE_TASK_ID
#!~~QUEUE TEMPLATE~~!#
# add method for setting TASK_ID and JOB_ID values for your queueing system
#elif [ "$USE_MY_QUEUE" = "true" ]; then
#	JOB_ID=$(get_my_queue_job_id)
#	TASK_ID=$(get_my_queue_task_id)
fi

# initialize all our important variables & directories
JOB_DIR=${SHRTCACHE}/$(whoami)/${JOB_ID}_${TASK_ID}
if [ "$USE_DB2_TGZ" = "true" ]; then
	batchsize=${USE_DB2_TGZ_BATCH_SIZE-1}
elif [ "$USE_DB2" = "true" ]; then
	batchsize=${USE_DB2_BATCH_SIZE-100}
fi
offset=$((batchsize*TASK_ID))
echo $offset $batchsize
INPUT_FILES=$(head -n $offset $EXPORT_DEST/file_list | tail -n $batchsize)

# log information about this job
log host=$(hostname)
log user=$(whoami)
log EXPORT_DEST=$EXPORT_DEST
log INPUT_SOURCE=$INPUT_SOURCE
log DOCKFILES=$DOCKFILES
log DOCKEXEC=$DOCKEXEC
log SHRTCACHE=$SHRTCACHE
log LONGCACHE=$LONGCACHE
log JOB_ID=$JOB_ID
log SGE_TASK_ID=$SGE_TASK_ID
log SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID
log TASK_ID=$TASK_ID

# validate required environmental variables
first=
fail=
for var in EXPORT_DEST INPUT_SOURCE DOCKFILES DOCKEXEC SHRTCACHE LONGCACHE JOB_ID TASK_ID; do
  if [ -z ${!var} ]; then
    if [ -z $first ]; then
      echo "the following required parameters are not defined: "
      fail=t
      first=f
    fi
    echo $var
  fi
  ! [ -z $fail ] && exit 1
done

JOB_DIR=${SHRTCACHE}/$(whoami)/${JOB_ID}_${TASK_ID}

if ! [ "$USE_DB2" = "false" ]; then
	OUTPUT=$EXPORT_DEST/$TASK_ID
else
	OUTPUT=${EXPORT_DEST}/$(sed "${TASK_ID}q;d" $EXPORT_DEST/joblist | awk '{print $2}')
fi
log OUTPUT=$OUTPUT
log INPUT_FILES=$INPUT_FILES
log JOB_DIR=$JOB_DIR

LOG_OUT=${SHRTCACHE}/rundock_${JOB_ID}_${SGE_TASK_ID}.out
LOG_ERR=${SHRTCACHE}/rundock_${JOB_ID}_${SGE_TASK_ID}.err

# bring directories into existence
mkdir -p $JOB_DIR/working
mkdir -p $DOCKFILES_TEMP

mkdir -p $OUTPUT
chmod -R 777 $OUTPUT

#cp -a $DOCKFILES/. $DOCKFILES_TEMP

mkdir $JOB_DIR/dockfiles
pushd $DOCKFILES
for f in $(find .); do
	fp=$PWD/$f
	jp=$JOB_DIR/dockfiles/$f
	mkdir -p $(dirname $jp)
	ln -s $fp $jp
done
popd
rm $JOB_DIR/dockfiles/INDOCK

# import restart marker, if it exists
# tells this script to ignore SIGUSR1 interrupts
trap '' SIGUSR1

if [ -f $OUTPUT/restart ]; then
	cp $OUTPUT/restart $JOB_DIR/working/restart
fi

# cleanup will:
# 1. move results/restart marker to $OUTPUT (if no restart marker, remove it from $OUTPUT if present)
# 2. remove the working directory
function cleanup {
        nout=$(ls $OUTPUT | grep OUTDOCK | wc -l)

        if [ $nout -ne 0 ] && ! [ -f $OUTPUT/restart ]; then
                log "Something seems wrong, my output is already full but has no restart marker. Removing items present in output and replacing with my results."
                rm $OUTPUT/*
                nout=0
        fi

        cp -p $JOB_DIR/working/OUTDOCK $OUTPUT/OUTDOCK.$nout
        cp -p $JOB_DIR/working/test.mol2.gz $OUTPUT/test.mol2.gz.$nout

        if [ -f $JOB_DIR/working/restart ]; then
                mv $JOB_DIR/working/restart $OUTPUT/restart
        elif [ -f $OUTPUT/restart ]; then
                rm $OUTPUT/restart
        fi

        rm -rf $JOB_DIR
}

# on exit, clean up files etc.
# setting this trap is good, as otherwise an error might interrupt the cleanup or abort the program too early
trap cleanup EXIT

pushd $JOB_DIR
$DOCKEXEC # this will produce an OUTDOCK with the version number
vers="3.7"
if [ -z "$(head -n 1 OUTDOCK | grep 3.7)" ]; then
	vers="3.8"
fi
rm OUTDOCK
popd

FIXINDOCK_SCRIPT=$JOB_DIR/fixindock.py
printf "
#!/bin/python3
import sys, os

def format_indock_arg(label, value):
	return '{:30s}{}'.format(label, value) + chr(10)
def fix_dockfiles_path(path):
	path = path[path.find('dockfiles'):]
	return '../' + path

with open(sys.argv[1], 'r') as indock, open(sys.argv[2], 'w') as savedest:
	for i, line in enumerate(indock):
		if '_file' in line:
			tokens = line.strip().split()
			label, value = tokens[0], tokens[1]
			if label == 'ligand_atom_file':
				savedest.write(format_indock_arg(label, '-'))
			elif label == 'output_file_prefix':
				savedest.write(format_indock_arg(label, 'test.'))
			else:
				value = fix_dockfiles_path(value)
				savedest.write(format_indock_arg(label, value))
		else:
			if i == 0:
				savedest.write('DOCK $vers parameter')
				continue
			savedest.write(line)" > $FIXINDOCK_SCRIPT

python3 $FIXINDOCK_SCRIPT $DOCKFILES/INDOCK $JOB_DIR/dockfiles/INDOCK

TARSTREAM_SCRIPT=$JOB_DIR/tarstream.py
printf "
#!/bin/python3
import tarfile, sys, time, os, gzip
for filename in sys.argv[1:]:
	tfile = tarfile.open(filename, 'r:gz')
	for t in tfile:
		f = tfile.extractfile(t)
		if not f or not (t.name.endswith('db2.gz') or t.name.endswith('db2')):
			continue
		data = f.read()
		if data[0:2] == bytes([31, 139]):
			data = gzip.decompress(data)
		sys.stdout.write(data.decode('utf-8'))
sys.stdout.close()" > $TARSTREAM_SCRIPT

pushd $JOB_DIR/working > /dev/null 2>&1

awk_cmd='{if (substr($0, 1, 1) == "S" && NF==8 && length==47) print substr($0, 1, 35); else print $0}'

log "starting DOCK vers=$vers"
(
	if [ $USE_DB2_TGZ = "true" ]; then
		if [ $vers = "3.7" ]; then
			python3 $TARSTREAM_SCRIPT $INPUT_FILES | awk "$awk_cmd"
		else
			python3 $TARSTREAM_SCRIPT $INPUT_FILES
		fi
	elif [ $USE_DB2 = "true" ]; then
		if [ $vers = "3.7" ]; then
			zcat -f $INPUT_FILES | awk "$awk_cmd"
		else
			zcat -f $INPUT_FILES
		fi
	fi
) | /usr/bin/time -v -o $OUTPUT/perfstats $DOCKEXEC $JOB_DIR/dockfiles/INDOCK &
dockpid=$!

function notify_dock {
	echo "notifying dock!"
	kill -10 $dockpid
}

trap notify_dock SIGUSR1

wait $dockpid
sleep 5 # bash script seems to jump the gun and start cleanup prematurely when DOCK is interrupted. This is stupid but effective at preventing this

# don't feel like editing DOCK src to change the exit code generated on interrupt, instead grep OUTDOCK for the telltale message
sigusr1=`tail OUTDOCK | grep "interrupt signal detected since last ligand- initiating clean exit & save" | wc -l`

log "finished!"

popd > /dev/null 2>&1

if [ $sigusr1 -ne 0 ]; then
	echo "s_rt limit reached!"
fi
exit 0
