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
