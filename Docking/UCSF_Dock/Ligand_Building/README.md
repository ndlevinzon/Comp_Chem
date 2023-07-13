# The 3D Ligand Building Pipeline
The 3D pipeline is a collection of scripts and software packages that enable the massively parallel creation of dockable 3D molecules.

## EZ Setup

### BKS Cluster
```
source /nfs/soft/dock/versions/dock38/pipeline_3D_ligands/env.(sh|csh)
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

For example, if one of your nodes went down and caused a bunch of jobs to fail, it would be safe to re-run ./submit-all-jobs.bash to re-submit those jobs. (assuming there are no jobs for that file currently queued/running)
