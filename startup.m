%% Search paths
% TODO: adapt for cases where there is another matlab folder already in use
userpath( fileparts( mfilename('fullpath') ) );    
addpath( genpath( userpath ) );
rmpath( genpath(  [userpath '/.git'] ) );
%rmpath( genpath(  [userpath '/test'] ) );
rmpath( genpath(  [userpath '/old'] ) );

user = getenv('USER');
hostname = getenv('HOSTNAME');
fprintf( 'USER : %s', user );
fprintf( '\nHOSTNAME : %s', hostname );
fprintf( '\nuserpath : %s', userpath );

%% CUDA
fprintf( '\nCUDA_PATH : %s', getenv( 'CUDA_PATH' ) );

%% ASTRA
% ASTRA 1.9 local, compiled 2019-09-09
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

%% MATLAB path
fprintf( '\nMATLAB_PATH : %s', getenv( 'MATLAB_PATH' ) );

%% Git repository version
fprintf( '\nGit commit ID : %s', git_commit_id );

%% Set default color map to grayscale instead of jet
set(groot, 'DefaultFigureColormap', gray)
close all;

fprintf( '\n' )
