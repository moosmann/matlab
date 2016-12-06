%% Scripts and functions
[~, hostname] = unix('echo $HOSTNAME');
if strcmp( hostname(1:3), 'max')
    local_matlab_path = '/asap3/petra3/gpfs/common/p05/jm/matlab';
else
    local_matlab_path = [getenv('HOME') '/matlab/'];
end
addpath( genpath( local_matlab_path ) );
rmpath( genpath(  [local_matlab_path '.git'] ) );
rmpath( genpath(  [local_matlab_path 'test'] ) );
rmpath( genpath(  [local_matlab_path 'old'] ) );

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
