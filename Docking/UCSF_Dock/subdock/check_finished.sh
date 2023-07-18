#!/bin/bash

echo "Checking for failed and finished jobs..."

rm FINISHED NOT-FINISHED
for i in `cat dirlist`
do
if [ -d ${i} ]
then
#echo ${i}
        if [ -f ${i}/OUTDOCK.0 ]
        then
                if [ -f ${i}/OUTDOCK.all ]
                        then
                        rm ${i}/OUTDOCK.all
                        fi
                cat ${i}/OUTDOCK.* > ${i}/OUTDOCK.all
                a=`grep "we reached the end" ${i}/OUTDOCK.all | wc -l | awk '{print $1}'`
                if [ $a -ge 1 ]
                        then
                        echo ${i} >> FINISHED
                        else
                        echo ${i} >> NOT-FINISHED
                        fi
        else
        echo ${i} >> NOT-FINISHED
        fi
else
echo ${i} >> NOT-FINISHED
fi
done

echo "You have this many failed jobs:"
wc -l NOT-FINISHED

echo "Making a file for job errors in case needed..."

rm LIST_OF_ERRORS
while read line; do
   echo $line >> LIST_OF_ERRORS
   cat $line/*err >> LIST_OF_ERRORS
   echo " " >> LIST_OF_ERRORS
done<NOT-FINISHED
echo "Errors for all failed jobs located in LIST_OF_ERRORS"

