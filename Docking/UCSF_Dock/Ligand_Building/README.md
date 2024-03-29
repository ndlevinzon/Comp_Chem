# The 3D Ligand Building Pipeline
The 3D pipeline is a collection of scripts and software packages that enable the massively parallel creation of dockable 3D molecules.

## Installing the 3D Pipeline
Installation of the 3D pipeline is somewhat tricky- the environment is very particular about which versions of which software should be used. It is thus easiest to copy the exact software packages we use from our servers- this should work provided they are installed in a 64-bit linux architecture, though depending on the distribution certain shared libraries may be missing.

### Setting up the installation root

First find a suitable directory which will serve as the root of the installation- this should be a directory that is visible from all nodes in the cluster.

Within this directory, create two sub-directories:
```
$ROOT_DIR/
    soft
    licenses
```
Copy or link your openeye and chemaxon licenses to the licenses directory- name them ".oe-license.txt" and ".jchem-license.cxl", respectively.

Next, clone the submission scripts from github to this directory.
```
git clone https://github.com/docking-org/zinc22-3d-submit
```

I like to rename this repository directory to just "submit", leaving the installation looking like this:
```
$ROOT_DIR/
    soft
    licenses/
        .oe-license.txt
        .jchem-license.cxl
    submit
```

Finally, create the "env.sh" and "env.csh" scripts in the ROOT_DIR as follows:

* bash
```
#!/bin/bash
# env.sh
base=<<ROOT_DIR>>
export BINDIR=$base/submit
export SHRTCACHE=<<TEMPDIR 1>>
export LONGCACHE=<<TEMPDIR 2>>
export SOFT_HOME=$base/soft
export LICENSE_HOME=$base/licenses
export PATH=$PATH:$base/submit
```
* csh
```
#!/usr/bin/csh
# env.csh
set base=<<ROOT_DIR>>
setenv SHRTCACHE <<TEMPDIR 1>>
setenv LONGCACHE <<TEMPDIR 2>>
setenv BINDIR $base/submit
setenv SOFT_HOME $base/soft
setenv LICENSE_HOME $base/licenses
setenv PATH $PATH\:$base/submit
```

<<ROOT_DIR>> should be the installation directory you chose.

<<TEMPDIR 1>> and <<TEMPDIR 2>> should be temporary directories available to all nodes on the cluster- often in distributed computing environments there are special directories set aside for this purpose, e.g /scratch

<<TEMPDIR 1>> will be used for short-term storage of small job files- it is thus appropriate to set this to a faster-access lower-capacity location, like /dev/shm. using /dev/shm can introduce problems, so it is safest to use the same value as LONGCACHE

<<TEMPDIR 2>> will be used for long-term storage of software files- thus it is more appropriate to set this to a slower-access higher-capacity location, like /tmp or /scratch.

## EZ Setup

### BKS Cluster
```
source /nfs/soft/dock/versions/dock38/pipeline_3D_ligands/env.(sh|csh)
```
* bash
```
#!/bin/bash
base=/nfs/soft/dock/versions/dock38/pipeline_3D_ligands
export BINDIR=$base/submit
export SHRTCACHE=/scratch
export LONGCACHE=/scratch
export SOFT_HOME=$base/soft
export LICENSE_HOME=$base/licenses
export PATH=$PATH:$base/submit
```
* csh
```
#!/usr/bin/csh

setenv SHRTCACHE /scratch
setenv LONGCACHE /scratch

set base=/nfs/soft/dock/versions/dock38/pipeline_3D_ligands

setenv BINDIR $base/submit
setenv SOFT_HOME $base/soft
setenv LICENSE_HOME $base/licenses
setenv PATH $PATH\:$base/submit
```

This environment will set up most of the required variables for you, as well as adds the submission scripts to your PATH, which means submission can be as simple as:

* bash
```
export INPUT_FILE=$HOME/myligands.smi
export OUTPUT_DEST=$HOME/myoutput
submit-all-jobs-slurm.bash
```
* csh
```
setenv INPUT_FILE $HOME/myligands.smi
setenv OUTPUT_DEST $HOME/myoutput
submit-all-jobs-slurm.bash
```

### Wynton Cluster
```
source /wynton/group/bks/soft/pipeline_3D_ligands/env.(sh|csh)
```

Similar to the BKS example, this environment will set most of the required variables for you

* bash
```
export INPUT_FILE=$HOME/myligands.smi
export OUTPUT_DEST=$HOME/myoutput
submit-all-jobs-sge.bash
```
* csh
```
setenv INPUT_FILE $HOME/myligands.smi
setenv OUTPUT_DEST $HOME/myoutput
submit-all-jobs-sge.bash
```

## Repackaging Output For Docking

The output of the 3D pipeline scripts will be a number of tar.gz files with roughly LINES_PER_JOB molecules contained per package. It is standard practice to repackage these smaller packages into larger packages for docking, as 50 molecules do not take long to process with DOCK.

### The DB2 File Format

DB2 files encode & compress 3D molecule rotomers, sometimes referred to as conformations. DB2 files do not store the dihedral angles of rotomers, but instead the per-atom coordinates resulting from those angles. Atoms of different rotomers whose coordinates overlap are merged together in db2 files (how this is done will be detailed more later), thus allowing DOCK to avoid double-calculating per-atom energies.

DB2 is a text format, with each row of information being printed on a new line (not to exceed 80 characters per line, due to old fortran nonsense). Mutiple DB2 entries can be stored in a db2 file- there is no limit on the number of entries that can be fit in a single db2 file. Thus concatenating any number of db2 files together still results in a valid db2 file.

There are seven(/eight) species of lines in a db2 entry- they will be listed in the order they appear.

* "M" lines serve a few functions. First is to define information about a db2 molecule entry- molecule name, smiles, properties, etc. The second is to create a boundary between distinct db2 entries

* "A" and "B" lines define atoms and bonds for a molecule. These are practically identical to atom/bond lines in .mol2 files- not much else to say here

* "X" lines define the set of all coordinates for the entry- these will be referenced later on. coordinates are guaranteed to be "distinct", in the sense that no coordinate is within a threshold distance (typically 0.001A) of any other coordinate in the set of all coordinates

* "R" lines define the coordinates of the rigid component. this is the basis component of the molecule that does not move between rotomers, usually a benzene ring or the like. the choice of which structure should be the rigid component is arbitrary, and in fact usually a separate db2 entry is created for each possible rigid component (usually each ring system) in a molecule

* "C" lines are "confs", not to be confused with conformations. We will get back to these.

* "S" lines or "set" lines define the rotomers/conformations. Each set entry may span multiple "S" lines to fit line-width requirements, and each set entry refers to a single rotomer/conformation. We will get back to these as well.

* "D" lines are "clusters", and are a feature of older db2 files. The cluster lines were meant to group together similar rotomers, though the idea was scrapped and newer db2 files no longer have them, thus "D" may now stand for "deprecated"
### What is the deal with sets and confs?

First off- the naming is super confusing. Confs don't describe conformations(rotomers), rather they describe a subset of a conformation. "Set" is an annoyingly general term for something that is so specific. These are the terms we are stuck with though! If I could turn back time, I would probably rename set->conformation and conf->subconformation

Sets are collections of confs, and confs are collections of atom xyz positions. Each set describes a conformation(rotomer) of the molecule. A single conf may be a part of one or more sets, thus confs can be used to define exactly where sets overlap. It is common to see confs comprised of a single atom- this may describe a situation where two rotomers happen to overlap at a particular atom. A conf may also be a part of just one set- this describes a situation where none of the atoms in the conf overlap with other sets, and is quite common.

### Example

Learning by example is best, so let's imagine the following (fake) rotomer set for a very simple molecule
```
+-------+
|   H   |
|   |   |
|O==C   |
+-------+
|       |
|       |
|O==C--H|
+-------+
```

Let us say the C atom lies at (0, 0). A somewhat abbreviated db2 for the set would look like this:
```
M fake
A 1 O
A 2 C
A 3 H
B 1 1 2 2
B 2 2 3 1
X 1 1 1 -1.000 +0.000
X 2 2 1 +0.000 +0.000
X 3 3 2 +0.000 +1.000
X 4 3 3 +1.000 +0.000
R 1 7 -1.000 +0.000
R 2 7 +0.000 +0.000
C 1 1 2
C 2 3 3
C 3 4 4
S 1 1 2 0 0 +0.000 +0.000
S 1 2 1 2
S 2 1 2 0 0 +0.000 +0.000
S 1 2 1 3
```
### Repackaging DB2 DOCK38
The following is a script for repackaging 3D pipeline results. First, here is the script:

```
#!/bin/bash
# make_tarballs.bash

# required parameter
TARBALL_SOURCE=$1
TARBALL_REPACK_DEST=$2

TARBALL_SOURCE=$(realpath $TARBALL_SOURCE)
TARBALL_REPACK_DEST=$(realpath $TARBALL_REPACK_DEST)

[ -z $TARBALL_SOURCE ] && echo "need to provide TARBALL_SOURCE as 1st arg!" && exit 1
[ -z $TARBALL_REPACK_DEST ] && echo "need to provide TARBALL_REPACK_DEST as 2nd arg!" && exit 1

# optional parameters
WORKING_DIRECTORY=${WORKING_DIRECTORY-/tmp/$(whoami)}
PACKAGES_PER_PACKAGE=${PACKAGES_PER_PACKAGE-100}
PACKAGE_TYPE=${PACKAGE_TYPE-db2.gz}
PACKAGE_TYPE_SHORT=$(echo $PACKAGE_TYPE | cut -d'.' -f1)

echo WORKING_DIRECTORY=$WORKING_DIRECTORY
mkdir -p $WORKING_DIRECTORY && cd $WORKING_DIRECTORY
mkdir -p output working tarball_split_list

echo finding
find $TARBALL_SOURCE -name '*.tar.gz' > tarball_list.txt
echo splitting
split -l $PACKAGES_PER_PACKAGE tarball_list.txt tarball_split_list/
echo working
cd working
for f in ../tarball_split_list/*; do
        for tb in $(cat $f); do
                ! [ -z $VERBOSE ] && echo tar --transform='s/^.*\///' -xf $tb '*.'$PACKAGE_TYPE 2>/dev/null
                tar --transform='s/^.*\///' -xf $tb '*.'$PACKAGE_TYPE 2>/dev/null
        done
        ! [ -z $VERBOSE ] && echo tar -czf $(basename $f).$PACKAGE_TYPE.tar.gz '*.'$PACKAGE_TYPE
        tar -czf $(basename $f).$PACKAGE_TYPE_SHORT.tar.gz *.$PACKAGE_TYPE
        mv $(basename $f).$PACKAGE_TYPE_SHORT.tar.gz $TARBALL_REPACK_DEST
        rm *.$PACKAGE_TYPE
        echo $(basename $f)
done
cd ..
rm -r $WORKING_DIRECTORY
echo Done! Results in $TARBALL_REPACK_DEST
```

Now, an example usage:

```
[user@gimel5 ~] bash make_tarballs.bash $PWD/H17P200_H19P400.smi.batch-3d.d/out $PWD/tarballs_repacked/H17P200_H19P400
finding
splitting
working
aa
ab
ac
ad
ae
af
ag
ah
ai
aj
Done! Results in $PWD/tarballs_repacked/H17P200_H19P400
```

It should be noted that this script will be effective for fairly small batches of molecules, e.g on the range of millions, rather than billions of molecules. For docking from ligands built using our pipeline with default options, running this script unmodified is sufficient for creating appropriately sized packages for docking. You may wish to edit WORKING_DIRECTORY to /scratch or some other larger directory if running out of space on /tmp is a concern. The /tmp directory typically only holds around 50G of data, which may not be enough for some workloads or environments.

## Script Arguments

Main submission scripts are named submit-all-jobs-slurm.bash and submit-all-jobs-sge.bash. These scripts use environment variables as arguments instead of usual command line ones.

E.g, on bash you would pass one of these arguments like so:
```
export INPUT_FILE=$PWD/example.smi
```
or on csh:
```
setenv INPUT_FILE $PWD/example.smi
```
Prior to running the script.

## Required Arguments

### INPUT_FILE

The input .smi file to be built. This file should contain only two columns of data: (SMILES, NAME) with no header.

### OUTPUT_DEST

The base directory for output to be stored. The script will create a sub-directory here named $INPUT_FILE.batch-3d.d

Within this output directory there are 3 sub-directories:
<ol>
<li> in </li>
<li> log </li>
<li> out </li>
</ol>

In contains the input file split into fragments and sub-fragments. By default the script first splits the input file into batches of 50000, then splits those batches into sub-batches of 50. Each individual job works on one of these sub-batches. Each array batch job works on one of the batches of 50000. All of the other directories alongside 'in' share the same directory structure.

Log contains log messages from the jobs. If you are re-submitting a file, be aware that log messages from previous runs on this file will be overwritten.

Out contains tar.gz output from each job. The tarballs should contain a number of 3d molecule formats for each molecule in the batch, including 1 or more db2.gz files.

### SOFT_HOME

Where software tarballs for the pipeline are stored. Symbolic links should be maintained in this directory according to the rules described in the "software arguments" section of this page. If you're sourcing a premade environment, don't worry about setting this value.

### LICENSE_HOME

Where software licenses are stored. Currently our licensed software includes jchem and openeye, licenses must be named .jchem-license.cxl and .oe-license.txt respectively. if you're sourcing a premade environment, don't worry about setting this value.

## Script Arguments

### SHRTCACHE

The base working directory for the script. By default it is /scratch

### LONGCACHE

The base directory for persistent files that are shared between jobs to go (i.e where software is installed). By default it is /scratch.

### CORINA_MAX_CONFS

How many nitrogen flapping configurations of each protomer corina should generate. By default only one is generated.

## Omega Arguments

These parameters correspond to torsion driving parameters described in the omega manual: https://docs.eyesopen.com/applications/omega/omega/omega_opt_params.html#torsion-driving-parameters

If you'd like to know more about how these parameters function, cross reference with the manual page.

### OMEGA_MAX_CONFS

Maximum configurations OMEGA will generate, default 600.

### OMEGA_ENERGY_WINDOW

Torsion energy window, if set to zero OMEGA will use an alternative rotatable bond dependent window method instead. Default is 12

### OMEGA_TORLIB

Torsion library- can choose between GubaV21 or Original, default is Original.

### OMEGA_FF

https://docs.eyesopen.com/toolkits/cpp/oefftk/OEFFConstants/OEMMFFSheffieldFFType.html#OEFF::OEMMFFSheffieldFFType::MMFF94Smod

Default is MMFF94Smod.

### OMEGA_RMSD

Sets rmsd for clustering and filtering conformations. If zero, omega will use an alternative rotatable-bond dependent method instead. Default is 0.5

## Job Submission Arguments

### SUBMIT_MODE

Choose the job submission method, choose between SGE, SLURM, or TEST_LOCAL. This will be automatically set if you use the job controller's corresponding superscript, e.g submit-all-jobs-slurm.bash. TEST_LOCAL will bypass the job controller and run the first input chunk in your shell.

### LINES_PER_BATCH

How many lines of the source .smi file should be processed per array batch job, default is 50000.

### LINES_PER_JOB

How many lines of the batch .smi file should be processed per array task, default is 50.

### MAX_BATCHES

Each batch job will contain LINES_PER_BATCH/LINES_PER_JOB jobs, and there will be a maximum of MAX_BATCHES batches submitted at any given time. By default this value is 25, which corresponds to 25,000 queued jobs at any given time if there are 1000 jobs per batch.

The submit-all script will block until less than MAX_BATCHES job arrays are in the queue. TODO: block until less than MAX_BATCHES total jobs are running or in the queue.

### SBATCH_ARGS

Additional arguments for the sbatch command. It is recommended to set a --time limit, as build jobs will save progress & terminate if they are still running two minutes before the --time limit.

### QSUB_ARGS

Additional arguments for the qsub command. Similar to slurm, it is recommended to set a time limit, but you will need to manually specify both s_rt & h_rt. In the example, we set s_rt to be a minute and thirty seconds before h_rt. s_rt is the point where jobs will save progress and terminate, h_rt is when they will be forcibly terminated, even if they've not finished saving.

## Software Options

All software variables will be set automatically if there exists a symbolic link in $SOFT_HOME matching the software variable's name, for example:
```
dock-latest -> DOCK.3.8.4.3d.tar.gz
jchem-latest -> jchem-19.15_r1.tar.gz
pyenv-latest -> lig_build_py3-3.7.1.tar.gz</nowiki>
```
They may also bet set manually- value is expected to be a path to a tar.gz file.

We use the following software:

* DOCK_VERSION

* JCHEM_VERSION

* PYENV_VERSION

* CORINA_VERSION

* OPENBABEL_VERSION

* EXTRALIBS_VERSION
  Note on EXTRALIBS- Run the pipeline with an empty EXTRALIBS package (but all other software accounted for) and see which shared libraries come up as missing in the error log. Locate all missing libraries and toss them in EXTRALIBS, they will be added to LD_LIBRARY_PATH

* JAVA_VERSION

## Examples

### Minimal Example
```
export INPUT_FILE=$PWD/example.smi
export OUTPUT_DEST=$PWD
bash submit-all-jobs-slurm.bash
```

### BKS Example - limit time to 2 hours, change batch size variables. Slurm tasks should automatically save progress when reaching their time limit.
```
export INPUT_FILE=$PWD/example.smi
export OUTPUT_DEST=$PWD/ligand_building
export SBATCH_ARGS="--time=02:00:00"
export LINES_PER_BATCH=20000
export LINES_PER_JOB=25
export MAX_BATCHES=15
bash submit-all-jobs-slurm.bash
```

### Wynton Example - limit time to 30 minutes, but set a soft limit 1:30 prior to the hard limit - the interrupt generated by the soft limit will signal the job to save progress for any resubmissions and exit.
```
export INPUT_FILE=$PWD/example.smi
export OUTPUT_DEST=$PWD/ligand_building
export QSUB_ARGS="-l s_rt=00:28:30 -l h_rt=00:30:00 -r y"
bash submit-all-jobs-sge.bash
```
## Resubmission

If your jobs for building have finished (or timed out), and you want to continue process whatever has not been processed yet, just run submit-all-jobs-slurm/sge again (with same env arguments). The submit-all script will detect which entries haven't finished and resubmit them.

## Repatriation

At BKS, we currently store the tarred output of the pipeline @ /nfs/exb/zinc22/tarballs. Currently, we use the following command to repatriate output from other clusters to our cluster:
```
### migrate_output.bash

for output in $OUTPUT_DEST/*.batch-3d.d; do
        echo "starting rsync on $output to $MIGRATE_USER@files2.docking.org"
        sshpass -f $PW_FILE rsync -arv $output/out $MIGRATE_USER@files2.docking.org:/nfs/exb/zinc22/tarballs/$(basename $output).out
done
```

sshpass is optional here but preferable for convenience's sake. Since files2.docking.org is only visible within the UCSF network, any clusters outside will need to maintain a network tunnel when rsyncing.

## Errors

Sometimes an output tarball will have few or no entries within. Certain molecule types will fail to be built, and often these molecules get bunched together (i.e if the input file is sorted by SMILES). Additionally, a small percentage of all molecules may fail to be processed by corina or amsol. If neither of these explain what is causing your missing entries, check that tarball's corresponding log entry for more info.

## Additional Notes

It is safe to re-run the same file multiple times- the script takes care of making sure not to re-run any jobs that have already completed successfully prior. This is only the case if that file's corresponding batch-3d.d output directory has not been moved or deleted.

== Developer Example: Building your own db2.tgz files & Submitting jobs ==

What you need:
# One or more db2(.gz) files
# An nfs directory(s) to store:
## docking input/output
## dockfiles
## dock executable
# subdock & rundock scripts

Create a list of all the db2 files you want to run docking against. The example below is merely a suggestion, make the list in any way you please so long as each entry is a *full* path (relative to your working directory) to a db2 (or db2.gz) file.

 <nowiki>
find $MY_DB2_SOURCE -type f -name "*.db2*" > my_db2_list</nowiki>

Split this list into reasonably sized chunks, our standard is 5000 but you can make them as large or small as you like. Do be careful about making the chunks larger- 5000 db2s is already quite heavy.

 <nowiki>
>> split -a 3 --lines=5000 my_db2_list db2_chunk.

>> ls
db2_chunk_aaa
db2_chunk_aab
db2_chunk_aac
...
my_db2_list</nowiki>

Create a db2.tgz archive from each of these lists.

 <nowiki>
for db2_chunk in db2_chunk.*; do
    tar -czf $db2_chunk.db2.tgz --files-from $db2_chunk
done</nowiki>


If you already have premade db2.tgz files, for example from the zinc22 3D archive, start the tutorial here.

Create a list of every db2.tgz archive. Each job will evaluate one db2.tgz archive. Again, this example is a suggestion applicable only if you've been following the tutorial up to this point.

 <nowiki>
find . -type f -name "db2_chunk.*.db2.tgz" > job_input_list</nowiki>

Now all you need to do is specify your docking parameters and launch the jobs. Set INPUT_SOURCE to be the job_input_list created in the previous step.

SGE

 <nowiki>
export DOCKEXEC=<dock executable path>
export DOCKFILES=<dockfiles path>
export EXPORT_DEST=<output directory path>
# optional arguments for the job controller. Note that these arguments are examples and not the only configuration recommended
export QSUB_ARGS="-l s_rt=00:28:00 -l h_rt=00:30:00 -l mem_free=2G"

export INPUT_SOURCE=job_input_list

bash <scripts directory>/sge/subdock.bash</nowiki>

SLURM

 <nowiki>
export DOCKEXEC=<dock executable path>
export DOCKFILES=<dockfiles path>
export EXPORT_DEST=<output directory path>
export SBATCH_ARGS="--time=00:30:00 --mem-per-cpu=2G"

export INPUT_SOURCE=job_input_list

bash <scripts directory>/slurm/subdock.bash</nowiki>

=== Large Docking Jobs ===

If your list of db2.tgz files is very large you may want to further split it. Each db2.tgz file in the job_input_list represents a job submitted to the queue, and often there is a limit on how many jobs can be queued at once.

In order to avoid this problem, we will need an automatic solution to split up our job_input_list and submit batches of jobs only when there is space left in the queue.

For example, imagine we want to submit in batches of 10,000 and limit total jobs to 50,000 (with a package size of 5000 this is 250M molecules in the queue maximum, submitting 50M at a time). The example I have shows how you would do this in slurm. 

 <nowiki>
#!/bin/bash
### submit_all_slurm.bash

BINDIR=$(dirname $0)

BATCH_SIZE=10000
MAX_QUEUED=50000

# the script is more portable if we provide the various parameters as arguments instead of hard-coding
INPUT_LIST=$1
BASE_EXPORT_DEST=$2 # this is the directory where further subdirectories will be created that contain docking job results
export DOCKEXEC=$3
export DOCKFILES=$4

# we can use our EXPORT_DEST as staging grounds for our input
mkdir -p $BASE_EXPORT_DEST/input

split --lines=$BATCH_SIZE -a 3 -n $INPUT_LIST $BASE_EXPORT_DEST/input/job_input.

export SBATCH_ARGS="--time=00:30:00 --mem-per-cpu=2G -J dock"

for job_input in $BASE_EXPORT_DEST/input/job_input.*; do

    export INPUT_SOURCE=$job_input
    input_num=$(printf $job_input | cut -d'.' -f2) # get the suffix of the split filename

    export EXPORT_DEST=$BASE_EXPORT_DEST/$input_num

    # loop forever
    while [ -z ]; do

        # counts how many jobs in total are pending or running on this user
        njobs=$(squeue -u $(whoami) -h -t pending,running -r | wc -l)
        # if you want to instead set a limit on how many *dock* jobs are pending or running you would just run the command through a filter
        # njobs=$(squeue -u $(whoami) -h -t pending,running -r | grep "dock" | wc -l)

        if [ $njobs -lt $((MAX_QUEUED-BATCH_SIZE)) ]; then
            break
        fi

        sleep 10
    done

    # the slurm subdock and rundock scripts need to live next to this script in a directory named "slurm"
    bash $BINDIR/slurm/subdock.bash
done</nowiki>

This script will run until all jobs are submitted, so for very large jobs you may want to keep the process alive in a screen.

== Tip: Using Wynton's $TMPDIR ==

<b>Doing this with SHRTCACHE_USE_ENV</b>

Before you run subdock, simply export the SHRTCACHE_USE_ENV option.

 <nowiki>
export SHRTCACHE_USE_ENV=TMPDIR</nowiki>

This will cause the script to use the $TMPDIR variable for SHRTCACHE.

== Example: Running a lot of docking jobs ==

* see [[ZINC22:Current status]] for more info about where ZINC can be found.

* 1. set up sdi files
 mkdir sdi
 export sdi=sdi
 ls /wynton/group/bks/zinc-22/H19/H19P0??/*.db2.tgz > $sdi/h19p0.in
 ls /wynton/group/bks/zinc-22/H19/H19P1??/*.db2.tgz > $sdi/h19p1.in
 ls /wynton/group/bks/zinc-22/H19/H19P2??/*.db2.tgz > $sdi/h19p2.in
 ls /wynton/group/bks/zinc-22/H19/H19P3??/*.db2.tgz > $sdi/h19p3.in
 and so on

* 2. set up INDOCK and dockfiles. rename dockfiles to dockfiles.$indockhash. On some nodes, the shasum command is called by sha1sum. Ultimately, renaming the dockfiles to a unique dockfiles is key. 

Note: As of 3/19/2021, this is no longer necessary

 bash
 indockhash=$(cat INDOCK | shasum | awk '{print substr($1, 1, 12)}')

* 3. super script:

 <nowiki>
export DOCKBASE=/wynton/group/bks/work/jji/DOCK
export DOCKFILES=$WORKDIR/dockfiles.21751f1bb16b
export DOCKEXEC=$DOCKBASE/docking/DOCK/bin/dock64
#export SHRTCACHE=/dev/shm # default
export SHRTCACHE=/scratch
export LONGCACHE=/scratch
export QSUB_ARGS="-l s_rt=00:28:00 -l h_rt=00:30:00 -l mem_free=2G"

for i in  sdi/*.in  ; do
        export k=$(basename $i .in)
	echo k $k
	export INPUT_SOURCE=$PWD/$i
	export EXPORT_DEST=$PWD/output/$k
	$DOCKBASE/docking/submit/sge/subdock.bash
done
</nowiki>

# 3a. to run for first time
 sh super

# 4. how to restart (to make sure complete, iterate until complete)

 sh super

# 5. check which output is valid (and broken or incomplete output)

# 6. extract all blazing fast

# 7. extract mol2

more soon, under active development, Jan 28.

== Appendix: Docking mono-cations of ZINC22 with DOCK3.8 on Wynton ==
Added by Ying 3/10/2021

To use: copy and paste the code section into terminal. '''Note to change the path where labelled with ''CHANGE this'' '''

* '''set up the folder to run docking. '''
Path to my example: /wynton/home/shoichetlab/yingyang/work/5HT-5a/10_AL-dock/zinc22_3d_build_3-10-2021
  mkdir zinc22_3d_build_3-10-2021
  cd zinc22_3d_build_3-10-2021

* '''copy INDOCK into dockfiles folder, and transfer to the created folder'''
  cp INDOCK dockfiles
  scp -r INDOCK dockfiles dt2.wynton.ucsf.edu:/path_to_created_folder

* '''get sdi of monocations of already built ZINC22 (<= H26 heavy atom count)'''
Modify to your own need...
 <nowiki>
mkdir sdi

foreach i (`seq 4 1 26`)
  set hac = `printf "H%02d" $i `
  echo $i $hac
  
  touch sdi/${hac}.sdi
  # CHANGE this: to your need
  foreach tgz (`ls /wynton/group/bks/zinc-22*/${hac}/${hac}[PM]???/*-O*.db2.tgz`)
    ls $tgz
    echo $tgz >> sdi/${hac}.sdi
  end
end
</nowiki>

* '''rename the dockfiles directory'''

Note: As of 3/19/2021 this step is no longer necessary

  indockhash=$(cat INDOCK | sha1sum | awk '{print substr($1, 1, 12)}')
  mv dockfiles dockfiles.${indockhash}

* '''write and run the super_run.sh'''
 <nowiki>
cat <<EOF > super_run.sh
export DOCKBASE=/wynton/group/bks/soft/DOCK-3.8.0.1
export DOCKEXEC=\$DOCKBASE/docking/DOCK/bin/dock64

# CHANGE here: path to the previously renamed dockfiles.\${indockhash}
### Note: as of 3/19/2021 renaming your dockfiles is no longer necessary
export DOCKFILES=/wynton/group/bks/work/yingyang/5HT-5a/10_AL-dock/zinc22_3d_build_3-10-2021/dockfiles.${indockhash}
export SHRTCACHE=/scratch
export LONGCACHE=/scratch
export QSUB_ARGS="-l s_rt=00:28:00 -l h_rt=00:30:00 -l mem_free=2G"

for i in  sdi/*.sdi  ; do
    export k=\$(basename \$i .sdi)
    echo k \$k
    export INPUT_SOURCE=$PWD/\$i
    export EXPORT_DEST=$PWD/output/\$k
    \$DOCKBASE/docking/submit/sge/subdock.bash
done
EOF

bash super_run.sh
</nowiki>

* '''keep submitting the super_run script until all db2s have been docked. '''
After all docking jobs finish, check the output. If no weird error, we can use a while loop to restart.
 <nowiki>
while true
do
  export jobN=$(qstat | grep -c 'rundock')
  if [[ $jobN -gt 0 ]] 
  then
    sleep 60
  else 
    bash super_run.sh
  fi
done
</nowiki>
When no new job is going to be submitted, use Ctrl+c to exit the while loop.

* '''extract scores from output. '''
 <nowiki>
cat << EOF > qsub_extract.csh
#\$ -S /bin/csh
#\$ -cwd
#\$ -pe smp 1
#\$ -l mem_free=100G
#\$ -l scratch=100G
#\$ -l h_rt=50:00:00
#\$ -j yes
#\$ -o extract_all.out

hostname
date

setenv DOCKBASE /wynton/group/bks/soft/DOCK-3.8.0.1

setenv dir_in $PWD

if ! (-d \$TMPDIR ) then
    if (-d /scratch ) then
        setenv TMPDIR /scratch/\$USER
    else
        setenv TMPDIR /tmp/\$USER
    endif
    mkdir -p \$TMPDIR
endif
pushd \$TMPDIR

ls -d \${dir_in}/output/*/*/ > dirlist

python \$DOCKBASE/analysis/extract_all_blazing_fast.py \
dirlist extract_all.txt -30

mv extract_all.* \$dir_in

popd

echo '---job info---'
qstat -j \$JOB_ID
echo '---complete---'
EOF

qsub qsub_extract.csh
</nowiki>

Another way is to run the command from the login node (Not recommended since sorting utilizes large memory)
 ls -d output/*/*/ > dirlist
 python $DOCKBASE/analysis/extract_all_blazing_fast.py dirlist extract_all.txt -20

* '''get poses in parallel'''
 <nowiki>
set score_file = $PWD/extract_all.sort.uniq.txt
set score_name = ${score_file:t:r}
set fileprefix = 'tmp_'
set number_per_file = 5000

set workdir  = $PWD/${score_name}_poses
mkdir -p $workdir 
cd $workdir

split --lines=$number_per_file --suffix-length=4 \
-d $score_file ${fileprefix}

set num  = ` ls ${fileprefix}* | wc -l `
echo "Number of score files to process:" $num

cat << EOF > qsub_poses.csh
#\$ -S /bin/csh
#\$ -cwd
#\$ -j yes
#\$ -pe smp 1
#\$ -l mem_free=5G
#\$ -l scratch=20G
#\$ -l h_rt=25:00:00
#\$ -t 1-$num

hostname
date

setenv DOCKBASE /wynton/group/bks/soft/DOCK-3.8.0.1

set list = \` ls \$PWD/${fileprefix}* \` 
set MOL = "\${list[\$SGE_TASK_ID]}"
set name = \${MOL:t:r}

python2 $DOCKBASE/analysis/getposes_blazing_faster.py \
"" \${MOL} $number_per_file poses_\${name}.mol2 test.mol2.gz

EOF

qsub qsub_poses.csh
cd ../
 </nowiki>

* '''Post-processing...'''

For example, if one of your nodes went down and caused a bunch of jobs to fail, it would be safe to re-run ./submit-all-jobs.bash to re-submit those jobs. (assuming there are no jobs for that file currently queued/running)
