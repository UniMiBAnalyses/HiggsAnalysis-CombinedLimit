pushd `dirname $BASH_SOURCE` > /dev/null
PWDX=`pwd`
popd > /dev/null

export PATH=${PATH}:${PWDX}/exe:${PWDX}/scripts
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PWDX}/lib
export PYTHONPATH=${PYTHONPATH}:${PWDX}/lib/python:${PWDX}/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PWDX}:/cvmfs/sft.cern.ch/lcg/releases/LCG_85swan3/vdt/0.3.6/x86_64-slc6-gcc49-opt/lib 

