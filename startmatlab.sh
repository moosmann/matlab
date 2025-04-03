#!/bin/zsh

# CUDA 11.8 for ASTRA
export PATH=/usr/local/cuda-11.8/bin:/software/cuda/cuda-11.8/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda-11.8/lib64:/usr/local/cuda-11.8/targets/x86_64-linux/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/software/cuda/cuda-11.8/targets/x86_64-linux/lib:/software/cuda/cuda-11.8/lib64:$LD_LIBRARY_PATH
export INCLUDE=/usr/local/cuda-11.8/include:/software/cuda/cuda-11.8/include:$INCLUDE
export CPATH=/usr/local/cuda-11.8/include:/software/cuda/cuda-11.8/include:$CPATH
export LIBRARY_PATH=/usr/local/cuda-11.8/lib64:/software/cuda/cuda-11.8/lib64:$LIBRARY_PATH
export LD_RUN_PATH=/usr/local/cuda-11.8/lib64:/software/cuda/cuda-11.8/lib64:$LD_RUN_PATH
export CUDNN_INCLUDE_PATH=/software/cuda/cuda-11.8/include
export CUDA_PATH=/software/cuda/cuda-11.8

# MATLAB
export MATLAB_VERSION=R2024b 
export MATLAB_USER_PATH=$PWD
export PATH=/opt/matlab/$MATLAB_VERSION/bin:$PATH

# ASTRA, maxwell installation requires CUDA 11.8 and MATLAB R2024
export ASTRA_PATH=/software/jupyter/.conda/envs/astra-toolbox-2.1.0
export ASTRA_SAMPLES_PATH=/asap3/petra3/gpfs/common/p05/astra/astra-toolbox/samples/matlab
export MATLABPATH=$ASTRA_PATH/matlab/mex:$ASTRA_PATH/matlab/tools:$MATLABPATH
export PATH=$ASTRA_PATH/bin:$PATH
#export LD_PRELOAD=/usr/lib64/libstdc++.so.6
#export LD_PRELOAD=/opt/nvidia/nsight-systems/2024.4.2/host-linux-x64/libstdc++.so.6:$LD_PRELOAD
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/nvidia/nsight-systems/2024.4.2/host-linux-x64:/usr/lib64:/usr/lib64

# Git
export GIT_COMMIT_ID=$(git rev-parse HEAD)

echo -e 'Set environment variables:'
echo -e 'TZ:' $TZ
echo -e 'PATH:' $PATH
echo -e 'LD_LIBRARY_PATH:' $LD_LIBRARY_PATH
echo -e 'LD_PRELOAD' $LD_PRELOAD
echo -e 'ASTRA_PATH:' $ASTRA_PATH
echo -e 'PYTHONPATH:' $PYTHONPATH
echo -e 'MATLAB_PATH:' $MATLAB_PATH
echo -e 'MATLABPATH:' $MATLABPATH
echo -e 'GIT_COMMIT_ID:' $GIT_COMMIT_ID
echo -e 'Starting MATLAB.\n'

#ulimit -u 128000
ulimit -u 16000
touch  $HOME/.matlab/startml
matlab_$MATLAB_VERSION -nosoftwareopengl 
