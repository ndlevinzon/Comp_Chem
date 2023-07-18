echo "Finding all test.mol2.gz.* to check..."

while read line; do
   ls -d $line/test.mol2.gz.* >> list_mol2_checking
done<$1

echo "Finding corrupt test.mol2.gz.* files..."

while read line; do
   echo $line 
   zcat $line >& temp 
   tail -1 temp >> temp.1
done<list_mol2_checking
grep gzip temp.1 > $2
rm temp temp.1
