%% Scripts and functions
%% TODO: adapt for cases where there is another matlab folder already in use
userpath( fileparts( mfilename('fullpath') ) );    
addpath( genpath( userpath ) );
rmpath( genpath(  [userpath '/.git'] ) );
%rmpath( genpath(  [userpath '/test'] ) );
rmpath( genpath(  [userpath '/old'] ) );

hostname = getenv('HOSTNAME');
fprintf( 'HOSTNAME : %s', hostname );
fprintf( '\nuserpath : %s', userpath );

%% CUDA
fprintf( '\nCUDA_PATH : %s', getenv( 'CUDA_PATH' ) );

%% ASTRA
% ASTRA 1.8 local
path_to_astra = '/asap3/petra3/gpfs/common/p05/astra/1.8/matlab';
% ASTRA 1.7.1 global, probably broken
%path_to_astra = '/opt/xray/astra-toolbox/1.7.1/matlab';
if strcmp( hostname(1:8), 'max-hzgg')    
    addpath( genpath( path_to_astra ) );
elseif strcmp( hostname(1:8), 'max-p3ag')   
    addpath( genpath( path_to_astra ) );
elseif strcmp( hostname(1:8), 'max-nova')
    %path_to_astra = '/asap3/petra3/gpfs/common/p05/astra/1.8_old/matlab';
    addpath( genpath( path_to_astra ) );
else    
    addpath( genpath( path_to_astra ) );
    %addpath( genpath( '/usr/share/astra/matlab' ) );
end
fprintf( '\nASTRA path : %s', path_to_astra );

%% ImageJ / Fiji
%path_to_imagej = '/asap3/petra3/gpfs/common/p05/jm/imagej/ImageJ.app/scripts';
%addpath( genpath( path_to_imagej ) );
%fprintf( '\nImageJ path : %s', path_to_imagej );
path_to_fiji = '/asap3/petra3/gpfs/common/p05/jm/fiji/Fiji.app/scripts';
addpath( genpath( path_to_fiji ) );
fprintf( '\nFiji path : %s', path_to_fiji );

%% Set default color map to grayscale instead of jet
set(groot, 'DefaultFigureColormap', gray)
close all;

%% BUG FIXES

% Fix error: An unexpected error occurred during CUDA execution. The
% CUDA error was: cannot set while device is active in this process.
%a = gpuArray(1); 
%clear a;
% Fix above bug for multi GPU system
for gpu_device_count = 1:gpuDeviceCount
    gpu = gpuDevice(gpu_device_count);
    a = gpuArray(1);
    fprintf( '\nGPU device index %u : %s', gpu_device_count, gpu.Name )
    clear gpu a
end

% Fix error: “dlopen: cannot load any more object with static TLS”
%ones(10) * ones(10);
