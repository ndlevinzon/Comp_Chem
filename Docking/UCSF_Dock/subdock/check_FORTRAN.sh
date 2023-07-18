#!/bin/bash

echo "Checking for failed and finished jobs..."

results_dir_name=$1

rm  NOT-FINISHED_${results_dir_name}
for i in `cat dirlist_${results_dir_name}`
do
if [ -d ${i} ]
then
#echo ${i}
	b=`grep "Warning: ieee_inexact is signaling" ${i}/0.err | wc -l | awk '{print $1}'`
	if [ $b -ge 1 ]
		then
		echo ${i} >> NOT-FINISHED_${results_dir_name}
		fi

else
echo ${i} >> NOT-FINISHED_${results_dir_name}
fi
done

echo "You have this many failed jobs:"
wc -l NOT-FINISHED_${results_dir_name}

echo "Making a file for job errors in case needed..."

rm LIST_OF_ERRORS_${results_dir_name}
while read line; do
   echo $line >> LIST_OF_ERRORS_${results_dir_name}
   cat $line/*err >> LIST_OF_ERRORS_${results_dir_name}
   echo " " >> LIST_OF_ERRORS_${results_dir_name}
done<NOT-FINISHED_${results_dir_name}
echo "Errors for all failed jobs located in LIST_OF_ERRORS"

echo "Now removing output from failed jobs before relaunching..."

while read line; do
   echo $line
   rm $line/*
done<NOT-FINISHED_${results_dir_name}

