%% Paths
% TODO: adapt for cases where there is another matlab folder already in use

% User path
startup_file =  mfilename('fullpath');
userpath( fileparts( startup_file ) );
% Search path
addpath( genpath( [ userpath filesep 'matlab' ] ) );

% Print info
user = getenv('USER');
hostname = getenv('HOSTNAME');
fprintf( 'USER : %s', user );
fprintf( '\nHOSTNAME : %s', hostname );
fprintf( '\nstartup file : %s.m', startup_file )
fprintf( '\nuserpath : %s', userpath );
fprintf( '\nCUDA_PATH : %s', getenv( 'CUDA_PATH' ) );

%% ASTRA
ASTRA_PATH = getenv( 'ASTRA_PATH' );
astra_path = [ ASTRA_PATH '/matlab' ];
fprintf( '\nASTRA_PATH : %s', ASTRA_PATH );
addpath( genpath( astra_path ) );
fprintf( '\nAdd ASTRA path : %s', astra_path );
astra_samples_path = getenv( 'ASTRA_SAMPLES_PATH' );
fprintf( '\nAdd ASTRA samples path : %s', astra_samples_path );
addpath( genpath( astra_samples_path ) );

%% MATLAB path
fprintf( '\nMATLAB_PATH : %s', getenv( 'MATLAB_PATH' ) );

%% ImageJ
%imagej_matlab = '/asap3/petra3/gpfs/common/p05/jm/imagej/Fiji.app/scripts';
%addpath(genpath(imagej_matlab)); % Update for your ImageJ2 (or Fiji) installation as appropriate ImageJ;

%% Git repository version
%fprintf( '\nGit commit ID : %s', git_commit_id );

%% Default figure properties
% Set default color map to grayscale instead of jet
set(groot, 'DefaultFigureColormap', gray)
set(groot,'DefaultFigureGraphicsSmoothing','off')
set( groot, 'DefaultFigureRenderer', 'painter')

%% Core info
%[~, ulim] = unix( 'ulimit -u;' );
fprintf( '\nulimit -u : ' )
unix( 'ulimit -u;' );
fprintf( evalc('feature(''numcores'');') );
%fprintf( '\n' )

%% Time to start
d = dir('~/.matlab/startml');
%t0 = second(datetime(d.datenum,'ConvertFrom','datenum'));
t0 = datetime(d.date);
t1 = datetime('now');
dt = t1 - t0;
fprintf('Startup time: %s',dt)

fprintf( '\n' )
