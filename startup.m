%% Scripts and functions
hostname = getenv('HOSTNAME');
if strcmp( hostname(1:3), 'max')
    %% !! PATH TO USER MATLAB FILES: ADJUST UNLESS YOU USE ~/matlab !!
    userpath( '/asap3/petra3/gpfs/common/p05/jm/matlab' );    
else
    userpath( [getenv('HOME') '/matlab'] );
end
addpath( genpath( userpath ) );
rmpath( genpath(  [userpath '/.git'] ) );
rmpath( genpath(  [userpath '/test'] ) );
rmpath( genpath(  [userpath '/old'] ) );

fprintf( 'HOSTNAME: %s\n', hostname );
fprintf( 'userpath: %s', userpath );

%% ASTRA 
% ASTRA 1.6
%addpath( genpath( '/usr/share/astra/matlab' ) );
% ASTRA 1.7.1
addpath( genpath( '/opt/xray/astra-toolbox/1.7.1/matlab' ) );

%% Set default color map to grayscale instead of jet
set(groot, 'DefaultFigureColormap', gray)
close all;

%% BUG FIXES

% Fix error: An unexpected error occurred during CUDA execution. The
% CUDA error was: cannot set while device is active in this process.
a = gpuArray(1); 
clear a;

% Fix error: “dlopen: cannot load any more object with static TLS”
%ones(10) * ones(10);
