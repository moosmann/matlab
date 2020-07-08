%% Paths
% TODO: adapt for cases where there is another matlab folder already in use

% User path
userpath( fileparts( mfilename('fullpath') ) );
% Search path
addpath( genpath( [ userpath filesep 'matlab' ] ) );

% Print info
user = getenv('USER');
hostname = getenv('HOSTNAME');
fprintf( 'USER : %s', user );
fprintf( '\nHOSTNAME : %s', hostname );
fprintf( '\nuserpath : %s', userpath );
fprintf( '\nCUDA_PATH : %s', getenv( 'CUDA_PATH' ) );

%% ASTRA
% ASTRA 1.99(?) local, compiled 2020-07-06
ASTRA_PATH = getenv( 'ASTRA_PATH' );
astra_path = [ ASTRA_PATH '/matlab' ];
addpath( genpath( astra_path ) );
fprintf( '\nASTRA_PATH : %s', ASTRA_PATH );
fprintf( '\nAdd ASTRA path : %s', astra_path );

%% MATLAB path
fprintf( '\nMATLAB_PATH : %s', getenv( 'MATLAB_PATH' ) );

%% Git repository version
%fprintf( '\nGit commit ID : %s', git_commit_id );

%% Set default color map to grayscale instead of jet
set(groot, 'DefaultFigureColormap', gray)

fprintf( '\n' )
