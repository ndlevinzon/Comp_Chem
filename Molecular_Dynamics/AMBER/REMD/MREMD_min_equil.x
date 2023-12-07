#!/bin/bash

if [ -f groupfile ]; then
  rm groupfile
fi

nrep=$(wc temperatures.dat | awk '{print $1}')
echo "Number of replicas: $nrep"
count=0

# Loop over each topology-coordinate pair
for FILENAME in $(ls *.rst7); do
  TOP_COORD_PAIR="${FILENAME%.rst7}"  # Remove the file extension
  TOP_FILE="${TOP_COORD_PAIR}.top"
  
  if [ -f "$TOP_FILE" ]; then
    for TEMP in $(cat temperatures.dat); do
      let count+=1
      REP=$(printf "%03d" $count)
      echo "TEMPERATURE: $TEMP K ==> FILE: equilibrate.mdin.$REP for $TOP_COORD_PAIR"

      # Modify equilibrate.mdin for the current temperature
      sed "s/XXXXX/$TEMP/g" equilibrate.mdin > temp
      sed "s/RANDOM_NUMBER/$RANDOM/g" temp > equilibrate.mdin.$REP

      # Append information to groupfile
      echo "-O -rem 0 -i equilibrate.mdin.$REP -o equilibrate.mdout.$REP -c $TOP_COORD_PAIR.rst7 -r equilibrate.rst.$REP -x equilibrate.mdcrd.$REP -inf equilibrate.mdinfo.$REP -p $TOP_FILE" >> groupfile

      rm -f temp
    done
    echo "#" >> groupfile
  else
    echo "Error: $TOP_FILE not found for $FILENAME."
  fi
done

echo "Number of REPLICAS = $nrep"
echo "Done."
