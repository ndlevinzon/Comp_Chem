#!/bin/bash

function log {
	echo "[$(date +%X)]: $@"
}

if [ -z $SHRTCACHE_USE_ENV ]; then
	SHRTCACHE=${SHRTCACHE-/dev/shm}
else
	SHRTCACHE=${!SHRTCACHE_USE_ENV}
fi

if [ "$USE_PARALLEL" = "true" ]; then
	TASK_ID=$1
elif [ "$USE_SLURM" = "true" ]; then
	JOB_ID=${SLURM_ARRAY_JOB_ID}
	TASK_ID=${SLURM_ARRAY_TASK_ID}
elif [ "$USE_SGE" = "true" ]; then
	JOB_ID=$JOB_ID
	TASK_ID=$SGE_TASK_ID
elif [ "$USE_CHARITY" = "true" ]; then
	JOB_ID=something
	TASK_ID=1
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

# charity engine has special overrides
if [ $USE_CHARITY = "true" ]; then
	DOCKFILES=/local/input/dockfiles
	tar -C /local/input -xf /local/input/dockfiles.tgz
	DOCKEXEC=/bin/dock64
	INPUT_SOURCE=/local/input
	EXPORT_DEST=/local/output
	USE_DB2_TGZ="true"
	batchsize=1
	INPUT_FILES=$(find /local/input -name '*.db2.tgz')
else
	TASK_ID_ACT=$(head -n $TASK_ID $EXPORT_DEST/joblist.$RESUBMIT_COUNT | tail -n 1)
	offset=$((batchsize*TASK_ID_ACT))
	INPUT_FILES=$(head -n $offset $EXPORT_DEST/file_list | tail -n $batchsize)
fi

# log information about this job
log host=$(hostname)
log user=$(whoami)
log EXPORT_DEST=$EXPORT_DEST
log INPUT_SOURCE=$INPUT_SOURCE
log DOCKFILES=$DOCKFILES
log DOCKEXEC=$DOCKEXEC
log SHRTCACHE=$SHRTCACHE
log JOB_ID=$JOB_ID
log TASK_ID=$TASK_ID
log TASK_ID_ACT=$TASK_ID_ACT

# validate required environmental variables
first=
fail=
for var in EXPORT_DEST INPUT_SOURCE DOCKFILES DOCKEXEC SHRTCACHE JOB_ID TASK_ID RESUBMIT_COUNT; do
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

OUTPUT=${EXPORT_DEST}/$TASK_ID_ACT
log OUTPUT=$OUTPUT
log INPUT_FILES=$INPUT_FILES
log JOB_DIR=$JOB_DIR
		
# bring directories into existence
mkdir -p $JOB_DIR/working
mkdir -p $OUTPUT

# we need to create this symbolic link in a directory relative to where DOCK will run
# the reason: fortran won't accept paths longer than 80 chars
# thus, all paths to grid files etc should be formatted like ../dockfiles/* to sidestep this path limit
ln -s $DOCKFILES $JOB_DIR/dockfiles

# import restart marker, if it exists
if [ -f $OUTPUT/restart ]; then
	cp $OUTPUT/restart $JOB_DIR/working/restart
fi

# cleanup will:
# 1. determine if any errors have occurred during processing
# 2. if no errors, move results/restart marker to $OUTPUT (if no restart marker & no errors, remove it from $OUTPUT if present)
# 2. remove the working directory
function cleanup {

	# don't feel like editing DOCK src to change the exit code generated on interrupt, instead grep OUTDOCK for the telltale message
	sigusr1=`tail $JOB_DIR/working/OUTDOCK | grep "interrupt signal detected since last ligand" | wc -l`
	complet=`tail $JOB_DIR/working/OUTDOCK | grep "close the file" | wc -l`
	nullres=`tail $JOB_DIR/working/OUTDOCK | grep "total number of hierarchies" | awk '{print $5}'`

	success=true
	if [ "$nullres" = "0" ]; then
		log "detected null result! your files may not exist or there was an error reading them"
		rm $JOB_DIR/working/*
		OUTPUT_SUFFIX=_nullres
		success=false
	fi
	if [ "$complet" = "0" ] && ! [ "$sigusr1" -ne 0 ]; then
		log "detected incomplete result!"
		OUTPUT_SUFFIX=_incomplete
		success=false
	fi
	if [ "$sigusr1" -ne 0 ]; then
		log "detected interrupt signal was received in OUTDOCK"
	fi

	nout=$RESUBMIT_COUNT
	if ! [ -f $OUTPUT/OUTDOCK.0 ] && [ "$success" = true ]; then
		nout=0 # if we haven't had a successful run yet, name it "0" regardless of RESUBMIT_COUNT, bloody confusing I know
		# otherwise SUBDOCK isn't quite sure if we've started/completed the run w/o listing contents of each output directory
		# make note of this scruple here
		echo "OUTDOCK.0 is OUTDOCK.$RESUBMIT_COUNT" > $OUTPUT/note
		# I *would* make a symbolic link, but then analysis scripts might double count the files
	fi

        cp $JOB_DIR/working/OUTDOCK $OUTPUT/OUTDOCK$OUTPUT_SUFFIX.$nout
        cp $JOB_DIR/working/test.mol2.gz $OUTPUT/test.mol2.gz$OUTPUT_SUFFIX.$nout

        if [ -f $JOB_DIR/working/restart ]; then
                mv $JOB_DIR/working/restart $OUTPUT/restart
        elif [ -f $OUTPUT/restart ] && ! [ "$complet" = "0" ] && ! [ "$nullres" = "0" ]; then
                rm $OUTPUT/restart
        fi

        rm -r $JOB_DIR
}

# on exit, clean up files etc.
# setting this trap is good, as otherwise an error might interrupt the cleanup or abort the program too early
trap cleanup EXIT

pushd $JOB_DIR 2>/dev/null 1>&2
$DOCKEXEC 2>/dev/null 1>&2 # this will produce an OUTDOCK with the version number
vers="3.7"
if [ -z "$(head -n 1 OUTDOCK | grep 3.7)" ]; then
	vers="3.8"
fi
rm OUTDOCK
popd 2>/dev/null 1>&2

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
				savedest.write('DOCK $vers parameter{}'.format(chr(10)))
				continue
			savedest.write(line)" > $FIXINDOCK_SCRIPT

python3 $FIXINDOCK_SCRIPT $DOCKFILES/INDOCK $JOB_DIR/INDOCK

TARSTREAM_SCRIPT=$JOB_DIR/tarstream.py
printf "
#!/bin/python3
import tarfile, sys, time, os, gzip, signal
from urllib.request import urlopen
signal.signal(signal.SIGTERM, signal.SIG_IGN)
signal.signal(signal.SIGUSR1, signal.SIG_IGN)
for filename in sys.argv[1:]:
	if filename.startswith('http://') or filename.startswith('https://'):
		fdobj = urlopen(filename)
		#fdobj.tell = lambda : 0 # lol
	else:
		fdobj = open(filename, 'rb')
	try:
		with fdobj, tarfile.open(mode='r|gz', fileobj=fdobj) as tfile:
			for t in tfile:
				f = tfile.extractfile(t)
				if not f or not (t.name.endswith('db2.gz') or t.name.endswith('db2')):
					continue
				data = f.read()
				if data[0:2] == bytes([31, 139]):
					data = gzip.decompress(data)
				sys.stdout.write(data.decode('utf-8'))
	except BrokenPipeError:
		# probably means DOCK was interrupted by timeout or crashed
		# printing this error message just clutters up the log, so swallow it
		# it can be determined if DOCK was interrupted or crashed from other markers
		pass
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
) | env time -v -o $OUTPUT/perfstats $DOCKEXEC $JOB_DIR/INDOCK 2>/dev/null &
dockppid=$!
# find actual DOCK PID by grepping for our executable in ps output, as well as the returned PID (which should be a parent to the actual DOCK process)
dockpid=$(ps -ef | awk '{print $8 "\t" $2 "\t" $3}' | grep $dockppid | grep $DOCKEXEC | awk '{print $2}')

function notify_dock {
	log "time limit reached- notifying dock!"
	kill -USR1 $dockpid
}

trap notify_dock SIGUSR1

sleeptime=2
footgun=0
while sleep $sleeptime && [ -z "$(kill -0 $dockppid 2>&1)" ]; do
	# protect people from this issue with tempconf stuff here
	footgun=$(tail OUTDOCK | grep "Warning. tempconf" | wc -l)
	if [ $footgun -gt 0 ]; then
		log "footgun alert! tempconf message detected- you seem to be using a DOCK executable that isn't compatible with 3.8 ligands!"
		log "see here: https://wiki.docking.org/index.php?title=SUBDOCK_DOCK3.8#Mixing_DOCK_3.7_and_DOCK_3.8_-_known_problems"
		log "going to kill DOCK and exit"
		kill -9 $dockpid
	fi
	sleeptime=5 # wait 2 seconds at first such that if the footgun issue arises, we don't give the process too much time to write out to disk. then wait 5 for subsequent loops
done

# if test.mol2.gz was produced, I guess keep it around. get rid of OUTDOCK though
if [ $footgun -gt 0 ]; then
	rm OUTDOCK
	echo "this problematic OUTDOCK was removed to save disk space. https://wiki.docking.org/index.php?title=SUBDOCK_DOCK3.8#Mixing_DOCK_3.7_and_DOCK_3.8_-_known_problems" > OUTDOCK
	echo "dockexec=$DOCKEXEC" >> OUTDOCK
fi
wait $dockppid
sleep 5 # bash script seems to jump the gun and start cleanup prematurely when DOCK is interrupted. This is stupid but effective at preventing this

log "finished!"

popd > /dev/null 2>&1

exit 0
