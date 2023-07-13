#!/bin/bash

sdi_file=$1

WORKDIR=$PWD

export DOCKFILES=$WORKDIR/dockfiles
export DOCKEXEC=/nfs/soft/dock/versions/dock38/executables/dock38_nogist
export SHRTCACHE=/scratch
export LONGCACHE=/scratch

export TMPDIR=/scratch
export SBATCH_EXEC=/usr/bin/sbatch
export SQUEUE_EXEC=/usr/bin/squeue
export SBATCH_ARGS="--time=19:28:00"

export k=$(basename ${sdi_file} .sdi)
echo k $k
export INPUT_SOURCE=$WORKDIR/$sdi_file
export EXPORT_DEST=$WORKDIR/output_${k}/$k

sh /nfs/soft/dock/versions/dock38/DOCK/ucsfdock/docking/submit/slurm/subdock.bash
