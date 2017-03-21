%% Scripts and functions
%% TODO: adapt for cases where there is another matlab folder already in use
userpath( fileparts( mfilename('fullpath') ) );    
addpath( genpath( userpath ) );
rmpath( genpath(  [userpath '/.git'] ) );
%rmpath( genpath(  [userpath '/test'] ) );
rmpath( genpath(  [userpath '/old'] ) );

hostname = getenv('HOSTNAME');
fprintf( 'HOSTNAME: %s\n', hostname );
fprintf( 'userpath: %s', userpath );

%% ASTRA
% ASTRA 1.8 local
path_to_astra = '/asap3/petra3/gpfs/common/p05/astra/1.8/matlab';
% ASTRA 1.7.1 global, probably broken
%path_to_astra = '/opt/xray/astra-toolbox/1.7.1/matlab';
if strcmp( hostname(1:8), 'max-hzgg')    
    addpath( genpath( path_to_astra ) );
elseif strcmp( hostname(1:8), 'max-p3ag')   
    addpath( genpath( path_to_astra ) );
else    
    addpath( genpath( path_to_astra ) );
    %addpath( genpath( '/usr/share/astra/matlab' ) );
end

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
