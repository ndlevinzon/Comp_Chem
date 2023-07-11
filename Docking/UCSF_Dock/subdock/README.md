Important note- although DOCK 3.8 is in the header of this article, SUBDOCK is perfectly capable of running DOCK 3.7 workloads, though some features of DOCK 3.8 will not be taken advantage of.

# Installing
```
git clone https://github.com/docking-org/SUBDOCK.git</nowiki>
```
IMPORTANT: subdock.bash expects to live in the same directory as rundock.bash!!!

* subdock.bash is located @ subdock.bash relative to the repository root.

* subdock.bash can be called directly from any location- it is not sensitive to the current working directory.

# What's New?

Compared to older scripts, SUBDOCK is easier to use, has more features, and is much more flexible!

## December 2022

* All jobs platforms (e.g slurm, sge) are supported on the same script

* GNU Parallel is now supported as a jobs platform! Ideal for small-scale local testing. https://www.gnu.org/software/parallel/

* Subdock can now be run on both db2.gz individual files & db2.tgz packages. A batch_size can be set for both types, allowing for more flexibility.

* Arguments can be provided environmentally, e.g "export KEY=VALUE" or on the command line e.g "--key=value"

* Subdock now prints out a superscript to copy-paste on success, convenient for re-submission.

* Fully restartable on all jobs platforms! See below section for an explanation on what this means, why it matters, and instructions on usage.

* INDOCK version header is automatically corrected, as are any file paths referenced by INDOCK.

## May 2023

* You can provide http(s) URLs to dockable files as your input in lieu of file paths!

* Charity engine is now supported as a jobs platform! More instructions for using this further down. (https://www.charityengine.com/)

* Subdock will automatically detect if your jobs failed- no need to use an extra script to check if your jobs have actually finished or not

# Supported Platforms

There are four platforms currently supported:

* SLURM
* SGE (Sun Grid Engine)
** note for BKS lab: the SGE queue on gimel does not have python3, your jobs will not work!
* GNU Parallel (for local runs- ideal for testing)
* Charity Engine

One of these platforms must be specified- SLURM is the default. These platforms can be set by the

```
--use-slurm=true
--use-sge=true
--use-parallel=true
--use-charity=true
```
Arguments, respectively

## Using Charity Engine

To use charity engine, you must have access to an executable of the charity engine CLI, as well as GNU parallel.

Additionally, you must provide your charity authentication details in the form of the CHARITY_AUTHKEY or --charity-authkey variable.

WIP, more specific instructions to come.

# Supported File Types

DOCK can be run on individual db2.gz files or db2.tgz tar packages.

The file type can be specified via the --use-db2=true or --use-db2-tgz=true arguments. db2.tgz is the default

Each job dispatched by SUBDOCK will consume BATCH_SIZE files, where BATCH_SIZE is equal to --use-db2-batch-size or --use-db2-tgz-batch-size depending on which file type is chosen.

The number of jobs dispatched by SUBDOCK is equal to ceil(N / BATCH_SIZE), where N is the total number of input files.

# Restartability

'''ONLY APPLICABLE FOR DOCK 3.8+!'''

Restartability means that we can impose arbitrary time limits on how long our jobs can run *without* losing our progress. Time limits can be as large or as small as we want them to be, even as little as a few minutes per job! This flexibility lets docking jobs efficiently fill in the gaps between longer-running jobs on the same ecosystem, thus they will be preferentially treated by whichever system is in charge of scheduling.

#How to use for your Job Platform#

On SLURM, runtime can be defined with the "--time" argument, e.g:

```subdock.bash --use-slurm=true --use-slurm-args="--time=00:30:00"```

This will allow our job to run for 30 minutes before progress is saved & copied out.
On GNU parallel this is accomplished with "--timeout", e.g:

```subdock.bash --use-parallel=true --use-parallel-args="--timeout 1800"```

On SGE, the same can be achieved using the s_rt and h_rt parameters, e.g:

```subdock.bash --use-sge=true --use-sge-args="-l s_rt=00:29:30 -l h_rt=00:30:00"```
This tells SGE to warn the job 30 seconds prior to the 30 minute hard limit. 
GNU and SLURM platforms will provide a hard-coded 30 seconds notice, whereas this notice period must be manually defined for SGE jobs.

# How to continue jobs

Run subdock.bash again with the same parameters (particularly EXPORT_DEST, INPUT_SOURCE, USE_DB2, USE_DB2_TGZ, USE_DB2_BATCH_SIZE, and USE_DB2_TGZ_BATCH_SIZE) to restart your jobs! If you saved the superscript SUBDOCK spits out on successful submission, you can simply call that. 

You'll know there is no more work to be done if SUBDOCK prints "all N jobs complete!", SUBDOCK will also tell you what proportion of jobs have not yet completed on each submission.

Output files are appended with a suffix indicating how many times the docking task has been resubmitted, e.g OUTDOCK.0 for the first attempt, OUTDOCK.1 for the second, etc.

Be careful not to overlap your submissions- there are no guardrails in place to prevent this from happening if you are not careful.

# Full Example - All Steps

This example assumes you have access to a DOCK executable and an installed scheduling system (SGE/SLURM/Parallel), but nothing else.

1. Source subdock code from github
 
```
git clone https://github.com/docking-org/SUBDOCK.git</nowiki>
```

2. Fetch dockfiles from DUDE-Z- we will use DRD4 for this example. **note- SUBDOCK automatically detects your DOCK version & corrects the INDOCK header accordingly
```
wget -r --reject="index.html*" -nH --cut-dirs=2 -l1 --no-parent https://dudez.docking.org/DOCKING_GRIDS_AND_POSES/DRD4/dockfiles/</nowiki>
```
3a. Get db2 database subset sample via ZINC-22. Example provided below:
```
wget http://files.docking.org/zinc22/zinc-22l/H17/H17P050/a/H17P050-N-laa.db2.tgz
wget http://files.docking.org/zinc22/zinc-22l/H17/H17P050/a/H17P050-N-lab.db2.tgz
wget http://files.docking.org/zinc22/zinc-22l/H17/H17P050/a/H17P050-N-lac.db2.tgz</nowiki>
```

You can select a db2 database subset via cartblanche22.docking.org- for wget-able files, choose the DOCK37 (*.db2.tgz) format, with URL download type. Multiple download types are supported, for example if you are on Wynton you can download Wynton file paths- removing the need to download the files yourself.

3b. If you downloaded the db2.tgz files yourself, create an sdi.in file from your database subset, which will serve as a list of files to evaluate. For example:
```
find $PWD -type f -name '*.db2.tgz' > sdi.in
```

4. Export the parameters we just prepared as environment variables. ** You need a DOCK executable! This can be found via our download server if you have a license, otherwise lab members can directly pull https://github.com/docking-org/dock3.git. On BKS cluster, some curated executables have been prepared with labels @ /nfs/soft/dock/versions/dock38/executables. DOCK 3.7 executables may be found here as well!

```
export INPUT_SOURCE=$PWD/sdi.in
export EXPORT_DEST=$PWD/output
export DOCKFILES=$PWD/dockfiles
export DOCKEXEC=/nfs/soft/dock/versions/dock38/executables/dock38_nogist
```

5. Choose a platform. You must select only one platform - mixing and matching is not supported.
```
export USE_SLURM=true|...
export USE_SGE=true|...
export USE_PARALLEL=true|...
```

Any value other than exactly "true" will be interpreted as false.

6a. Run docking!
 <nowiki>
bash ~/SUBDOCK/subdock.bash</nowiki>

6b. You can also use command line arguments instead of environment export, if desired. These can be mixed and matched.
```
export DOCKEXEC=$PWD/DOCK/ucsfdock/docking/DOCK/dock64
bash ~/SUBDOCK/subdock.bash --input-source=$PWD/sdi.in --export-dest=$PWD/output --dockfiles=$PWD/dockfiles --use-slurm=true
```

7. After executing subdock, it will print out a convenient "superscript" to copy & paste, for any future re-submissions.

# Mixing DOCK 3.7 and DOCK 3.8 - known problems

'''Headline: Though SUBDOCK is compatible with DOCK 3.7, and will allow docking of ligands built for 3.8 in 3.7, it is NOT RECOMMENDED to do this without using a specially prepared 3.7 executable!'''

If you're running DOCK 3.8 against recently built ligands, you may encounter error messages that look like this:
 <nowiki>       1      2 bonds with error
Error. newlist is not big enough</nowiki>

Or worse, like this:
 <nowiki> Warning. tempconf = 0
         1597 ->            0 ->            0</nowiki>

The latter error messages have the potential to cause some serious damage, as they are emitted very frequently & may consume excessive disk space. SUBDOCK will check for these messages periodically during DOCK's runtime & kill the process if they are found.

If you are on 3.8 and are encountering these messages still, use the dock38_nogist executable described in [[How_to_install_DOCK_3.8#Prebuilt_Executable]]. This version voids the code related to the GIST scoring function, which is responsible for these errors.

If you are using 3.7 still, it is possible to prepare a version that keeps everything the same, except without the dangerous "tempconf" message.

== SUBDOCK help splash - all argument descriptions & defaults ==
 <nowiki>
[user@machine SUBDOCK]$ ./subdock.bash --help
SUBDOCK! Run docking workloads via job controller of your choice
# Required arguments
expected env arg: EXPORT_DEST, --export-dest
arg description: nfs output destination for OUTDOCK and test.mol2.gz files

expected env arg: INPUT_SOURCE, --input-source
arg description: nfs directory containing one or more .db2.tgz files OR a file containing a list of db2.tgz files

expected env arg: DOCKFILES, --dockfiles
arg description: nfs directory containing dock related files and INDOCK configuration for docking run

expected env arg: DOCKEXEC, --dockexec
arg description: nfs path to dock executable

# Job controller settings
optional env arg missing: USE_SLURM, --use-slurm
arg description: use slurm
defaulting to false

optional env arg missing: USE_SLURM_ARGS, --use-slurm-args
arg description: addtl arguments for SLURM sbatch command
defaulting to 

optional env arg missing: USE_SGE, --use-sge
arg description: use sge
defaulting to false

optional env arg missing: USE_SGE_ARGS, --use-sge-args
arg description: addtl arguments for SGE qsub command
defaulting to 

optional env arg missing: USE_PARALLEL, --use-parallel
arg description: use GNU parallel
defaulting to false

optional env arg missing: USE_PARALLEL_ARGS, --use-parallel-args
arg description: addtl arguments for GNU parallel command
defaulting to 

# input settings
optional env arg missing: USE_DB2_TGZ, --use-db2-tgz
arg description: dock db2.tgz tar files
defaulting to true

optional env arg missing: USE_DB2_TGZ_BATCH_SIZE, --use-db2-tgz-batch-size
arg description: how many db2.tgz to evaluate per batch
defaulting to 1

optional env arg missing: USE_DB2, --use-db2
arg description: dock db2.gz individual files
defaulting to false

optional env arg missing: USE_DB2_BATCH_SIZE, --use-db2-batch-size
arg description: how many db2.gz to evaluate per batch
defaulting to 100

# Addtl job configuration
optional env arg missing: MAX_PARALLEL, --max-parallel
arg description: max jobs allowed to run in parallel
defaulting to -1

optional env arg missing: SHRTCACHE, --shrtcache
arg description: temporary local storage for job files
defaulting to /scratch

optional env arg missing: LONGCACHE, --longcache
arg description: longer term storage for files shared between jobs
defaulting to /scratch
# Miscellaneous
optional env arg missing: SUBMIT_WAIT_TIME, --submit-wait-time
arg description: how many seconds to wait before submitting
defaulting to 5

optional env arg missing: USE_CACHED_SUBMIT_STATS, --use-cached-submit-stats
arg description: only check completion for jobs submitted in the latest iteration. Faster re-submission, but will ignore jobs that have been manually reset
defaulting to false
