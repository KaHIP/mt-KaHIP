#!/bin/bash

NCORES=4
unamestr=`uname`
if [[ "$unamestr" == "Linux" ]]; then
        NCORES=`grep -c ^processor /proc/cpuinfo`
fi

if [[ "$unamestr" == "Darwin" ]]; then
        NCORES=`sysctl -n hw.ncpu`
fi

rm -rf deploy
rm -rf build
mkdir build
cd build
cmake ../ -DCMAKE_BUILD_TYPE=Release $1
make -j $NCORES
cd ..

mkdir deploy
cp ./build/mtkahip deploy/
cp ./build/evaluator deploy/
cp ./build/graphchecker deploy/

cp ./build/libinterfacemtkahip_static.a deploy/libmtkahip.a
cp ./build/libinterfacemtkahip.so deploy/libmtkahip.so
cp ./interface/mtkahip_interface.h deploy/
