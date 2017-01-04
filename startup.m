%% Scripts and functions
%% TODO: adapt for cases where there is another matlab folder already in use
userpath( fileparts( mfilename('fullpath') ) );    
addpath( genpath( userpath ) );
rmpath( genpath(  [userpath '/.git'] ) );
rmpath( genpath(  [userpath '/test'] ) );
rmpath( genpath(  [userpath '/old'] ) );

fprintf( 'HOSTNAME: %s\n', getenv('HOSTNAME') );
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
