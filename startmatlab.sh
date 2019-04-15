#!/bin/sh

# Set timezone. If empty it can raise a Matlab warning when batch jops are executed
export TZ='Europe/Berlin'

# NVIDIA CUDA
export CUDA_PATH=/usr/local/cuda-9.1
export PYTHONPATH=$CUDA_PATH/lib64/:$PYTHONPATH

# MATLAB
export MATLAB_VERSION=R2017b

# ASTRA
export PATH=/opt/matlab/$MATLAB_VERSION/bin/:$PATH
export ASTRA_PATH=/asap3/petra3/gpfs/common/p05/astra/astra
export LD_LIBRARY_PATH=$ASTRA_PATH/lib:$CUDA_PATH/lib64/:$LD_LIBRARY_PATH
export MATLAB_PATH=$ASTRA_PATH/matlab/mex/:$ASTRA_PATH/matlab/tools/:$ASTRA_PATH/samples/matlab/:$MATLABPATH
export PYTHONPATH=$ASTRA_PATH/python:$PYTHONPATH

#export LD_LIBRARY_PATH=/usr/local/cuda/lib64/:$LD_LIBRARY_PATH

# local gcc 4.9, required to use matlab in python
#export LD_LIBRARY_PATH=/home/moosmanj/gcc/gcc/lib64/:$LD_LIBRARY_PATH
#export PATH=/home/moosmanj/gcc/gcc/bin:$PATH

# Git
export GIT_COMMIT_ID=$(git rev-parse HEAD)

echo -e 'Set environment variables:'
echo -e 'TZ:' $TZ
echo -e 'PATH:' $PATH
echo -e 'LD_LIBRARY_PATH:' $LD_LIBRARY_PATH
echo -e 'ASTRA_PATH:' $ASTRA_PATH
echo -e 'PYTHONPATH:' $PYTHONPATH
echo -e 'MATLAB_PATH:' $MATLAB_PATH
echo -e 'GIT_COMMIT_ID:' $GIT_COMMIT_ID
echo -e 'Starting MATLAB.\n'

matlab_$MATLAB_VERSION
