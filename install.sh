#/bin/sh

base=`pwd`

echo "start to install rbsSeeker..."

cd $base/thirdUtils/BamTools
make clean_api
make api
message=`make`

cd $base/thirdUtils/RNAfoldLib
make clean
message=`make`

cd $base/bioUtils
make clean
message=`make`

cd $base/src/rbsSeeker
make clean
message=`make`

cd $base

echo "installation finished!"
