#!/bin/bash

#  setup ROOT environment
source /usr/local/root/bin/thisroot.sh

#  drop_from_path taken from /root/build/bin/thisroot.sh
drop_from_path()
{
   # Assert that we got enough arguments
   if test $# -ne 2 ; then
      echo "drop_from_path: needs 2 arguments"
      return 1
   fi

   local p=$1
   local drop=$2

   newpath=`echo $p | sed -e "s;:${drop}:;:;g" \
                          -e "s;:${drop}\$;;g"   \
                          -e "s;^${drop}:;;g"   \
                          -e "s;^${drop}\$;;g"`
}


#  setup GEANT4 environment
cd /home/jdamgov/geant4.10.05.p01-install/bin/
/home/jdamgov/geant4.10.05.p01-install/bin/geant4.sh
cd -

export G4BASE=/home/jdamgov/geant4.10.05.p01/source
export G4INSTALL=/home/jdamgov/geant4.10.05.p01-install

source /home/jdamgov/geant4.10.05.p01-install/share/Geant4-10.5.1/geant4make/geant4make.sh
export G4BIN="$PWD"

export LD_LIBRARY_PATH=/home/jdamgov/geant4.10.05.p01-install/lib:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=/home/jdamgov/geant4.10.05.p01-install/lib:$DYLD_LIBRARY_PATH
export SHLIB_PATH=/home/jdamgov/geant4.10.05.p01-install/lib
export LIBPATH=/home/jdamgov/geant4.10.05.p01-install/lib

