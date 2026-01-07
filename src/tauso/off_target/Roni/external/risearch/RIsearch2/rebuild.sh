#!/bin/bash

if [ -f bin/risearch2.precompiled.x ];
then
   echo "bin/risearch2.precompiled.x exists."
else
   cp bin/risearch2.x bin/risearch2.precompiled.x
fi
rm -rf libdivsufsort-2.0.1/build || exit 1
mkdir -p libdivsufsort-2.0.1/build || exit 1
cd libdivsufsort-2.0.1/build || exit 1
cmake -DCMAKE_BUILD_TYPE="Release" -DBUILD_DIVSUFSORT64:BOOL=ON -DUSE_OPENMP:BOOL=ON -DBUILD_SHARED_LIBS:BOOL=OFF .. || exit 1
make || exit 1
cd - || exit 1
cd src || exit 1
make clean || exit 1
make || exit 1
cd -

echo
echo '-------------------------------------------------'
echo 'You should now be able to run RIsearch2 by typing'
echo '-------------------------------------------------'

echo "$PWD"/bin/risearch2.x


