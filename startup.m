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
% ASTRA 1.9 local
ASTRA_PATH = getenv( 'ASTRA_PATH' );
fprintf( '\nASTRA_PATH : %s', ASTRA_PATH );
astra_path = [ ASTRA_PATH '/matlab' ];
if strcmp( hostname(1:8), 'max-hzgg')    
    addpath( genpath( astra_path ) );
elseif strcmp( hostname(1:8), 'max-p3ag')   
    addpath( genpath( astra_path ) );
elseif strcmp( hostname(1:8), 'max-nova')
    addpath( genpath( astra_path ) );
else    
    addpath( genpath( astra_path ) );
end
fprintf( '\nAdd ASTRA path : %s', astra_path );

%% ImageJ / Fiji, currently not working well
%path_to_fiji = '/asap3/petra3/gpfs/common/p05/jm/fiji/Fiji.app/scripts';
%addpath( genpath( path_to_fiji ) );
%fprintf( '\nFiji path : %s', path_to_fiji );

%% MATLAB path
fprintf( '\nMATLAB_PATH : %s', getenv( 'MATLAB_PATH' ) );

%% Git repository version
fprintf( '\nGit commit ID : %s', git_commit_id );

%% Set default color map to grayscale instead of jet
set(groot, 'DefaultFigureColormap', gray)
close all;

%% BUG FIXES

% Fix error: An unexpected error occurred during CUDA execution. The
% CUDA error was: cannot set while device is active in this process.
%a = gpuArray(1); 
%clear a;
% Fix above bug for multi GPU system
% for gpu_device_count = 1:gpuDeviceCount
%     gpu = gpuDevice(gpu_device_count);
%     a = gpuArray(1);
%     fprintf( '\nGPU device index %u : %s', gpu_device_count, gpu.Name )
%     clear gpu a
% end

% Fix error: “dlopen: cannot load any more object with static TLS”
%ones(10) * ones(10);

fprintf( '\n' )
