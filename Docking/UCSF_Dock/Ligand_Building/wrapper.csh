#! /bin/csh

# source enviorment. 
source ~stefan/.cshrc_dbgen_corina
source /nfs/soft/jchem/jchem-19.15/env.csh
#set number_per_db2 = 1000
#set number_per_db2 = 5
set number_per_db2 = 10
set ph = 7.4

set file = $1
set fileprefix = ${file:r}


set pwd = `pwd`
set pathdir = $pwd
set workdir = $pwd/${fileprefix}

echo ${workdir}
echo ${file} ${fileprefix}

if (-e ${workdir}) then
 echo "${workdir} exists"
 exit
endif

mkdir ${workdir}
cd ${workdir}

mkdir ${workdir}/logs

ln -s ../${file} .

echo "split --lines=$number_per_db2 ${file}"

split --lines=$number_per_db2 ${file} --suffix-length=6


foreach splitfile ( ` ls x* | grep -v db2.gz ` ) 
echo ${splitfile}

# make sure that the link is pointing to something.  
#set lsoutput = `ls -l sgejob_*/${splitfile}.db2.gz`
#echo "WHAT:: $lsoutput"
#if ("$lsoutput" == "") then
#    rm ${splitfile}.db2.gz
#endif

if (-e ${splitfile}.db2.gz) then
#   echo "I AM HERE(2)"
   echo "${splitfile} has been submited for generations." 
   continue
endif

#rm stdout_${splitfile} stderr_${splitfile} script_qsub_${splitfile}.csh 
#echo "I AM HERE(3)"

cat << EOF >! qsub_${splitfile}.csh
#\$ -S /bin/csh
#\$ -cwd
#\$ -q all.q
#\$ -o ${workdir}/logs/stdout_${splitfile}
#\$ -e ${workdir}/logs/stderr_${splitfile}
#\$ -l h=!het&!n-1-28&!n-1-124&!gimel&!he

# source enviorment. 
source ~stefan/.cshrc_dbgen_corina
source /nfs/soft/jchem/jchem-19.15/env.csh
echo $PATH
hostname
date

set SCRATCH_DIR = /scratch
if ! (-d \$SCRATCH_DIR ) then
    SCRATCH_DIR=/scratch
endif
set username = `whoami`

set TASK_DIR = "\$SCRATCH_DIR/\${username}/\$JOB_ID"
echo \$TASK_DIR

mkdir -p \${TASK_DIR}
cd \${TASK_DIR}
pwd

cp ${workdir}/${splitfile} .

# note that the pining script's inputs should not be in quotes ('' or "").
/nfs/home/tbalius/zzz.github/DOCK/common/on-one-core - ${DOCKBASE}/ligand/generate/build_database_ligand.sh -H $ph ${splitfile} --no-db --save-table

cd ${workdir}
mkdir ${workdir}/${splitfile}_build

echo copying
ls -l \${TASK_DIR}/finished/*/*.db2.gz
ls -l \${TASK_DIR}/finished/*/*.db.gz 

# copy finshed directory
mv \${TASK_DIR}/finished/ ${workdir}/${splitfile}_build/.
rm -r \${TASK_DIR}

#mv ${splitfile} ${workdir}/${splitfile}_build/.

which cxcalc

EOF

set name = `whoami`

while ( `qstat -u ${name} | wc -l ` > 4950 )
  sleep 10
end

qsub qsub_${splitfile}.csh

#exit

end
