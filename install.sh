#!/bin/bash
# for use in bigdata machine only
rm -r build
mkdir build
cd build

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/local/opt/libomp/lib 
export LIBRARY_PATH=${LIBRARY_PATH}:/usr/local/opt/libomp/lib 
cmake ../source -DCMAKE_INSTALL_PREFIX=${LOCAL} -DCMAKE_INSTALL_LIBDIR=${LOCAL}/lib
make install