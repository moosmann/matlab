%% Scripts and functions
addpath( genpath( [getenv('HOME') '/matlab'] ) );
rmpath( genpath(  [getenv('HOME') '/matlab/.git'] ) );
rmpath( genpath(  [getenv('HOME') '/matlab/test'] ) );
rmpath( genpath(  [getenv('HOME') '/matlab/old'] ) );

%% ASTRA 
% ASTRA 1.6
%addpath( genpath( '/usr/share/astra/matlab' ) );
% ASTRA 1.7.1
addpath( genpath( '/opt/xray/astra-toolbox/1.7.1/matlab' ) );

%% Set default color map to grayscale instead of jet
set(groot, 'DefaultFigureColormap', gray)
close all;

%% BUG FIXES

% Shall fix error: An unexpected error occurred during CUDA execution. The
% CUDA error was: cannot set while device is active in this process.
a = gpuArray(1); 
clear a;

% Fix Matlab bug: “dlopen: cannot load any more object with static TLS”
%ones(10) * ones(10);
