# build all libraries and mex files
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug
make VERBOSE=0
cd ..
# copy mex files and libraries into bin directory
echo "copying mex files to final MEXs directory"
mkdir ../MEXs
rm -r ../MEXs/*
cp build/MEXs/*.mexa64 ../MEXs
# copy libraries needed by mex files
#echo "copying libs needed by mex files"
#cp build/libs/*.so ../libs
echo ""
