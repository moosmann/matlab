% TODO: adapt for cases where there is another matlab folder already in use

%% User search path
startup_file =  mfilename('fullpath');
userpath( fileparts( startup_file ) );
addpath( genpath( [ userpath filesep 'matlab'] ) );

%% ASTRA search path
ASTRA_PATH = getenv('ASTRA_PATH');
astra_path = [ ASTRA_PATH '/matlab'];
addpath( genpath( astra_path ) );
%astra_samples_path = getenv('ASTRA_SAMPLES_PATH');
%addpath( genpath( astra_samples_path ) );

if isempty(getCurrentTask())

    %% User info
    user = getenv('USER');
    hostname = getenv('HOSTNAME');
    fprintf('USER : %s', user );
    fprintf('\nHOSTNAME : %s', hostname );
    fprintf('\nstartup file : %s.m', startup_file )
    fprintf('\nuserpath : %s', userpath );
    fprintf('\nCUDA_PATH : %s', getenv('CUDA_PATH') );

    %% ASTRA info
    fprintf('\nASTRA_PATH : %s', ASTRA_PATH );
    fprintf('\nAdd ASTRA path : %s', astra_path );
 %   fprintf('\nAdd ASTRA samples path : %s', astra_samples_path );

    %% MATLAB path
    fprintf('\nMATLAB_PATH : %s', getenv('MATLAB_PATH') );
    fprintf('\nMATLABPATH : %s', getenv('MATLABPATH') );
    fprintf('\nMATLAB_USER_PATH : %s', getenv('MATLAB_USER_PATH') );

    %% ImageJ
    fprintf('\nIMAGEJ : %s', getenv('IMAGEJ') );
    fprintf('\nIMAGEJ_MACROS : %s', getenv('IMAGEJ_MACROS') );
    %addpath(genpath(imagej_matlab)); % Update for your ImageJ2 (or Fiji) installation as appropriate ImageJ;

    %% Default figure properties
    % Set default color map to grayscale instead of jet
    %set(groot,'DefaultFigureColormap', gray)
    %set(groot,'DefaultFigureColormap', stern_special)
    %set(groot,'DefaultFigureGraphicsSmoothing','off')
    %set( groot,'DefaultFigureRenderer','painter')

    %% Number of cores
    %[~, ulim] = unix('ulimit -u;');
    fprintf('\nulimit -u :')
    unix('ulimit -u;');
    fprintf( evalc('feature(''numcores'');') );
    %fprintf('\n')

    % %% Render info
    % fprintf('Rendererinfo:\n')
    % disp(rendererinfo)

        %% Git repository version
    %fprintf('\nGit commit ID : %s', git_commit_id );

    %% Time to start
    d = dir('~/.matlab/startml');
    t0 = datetime(d.date);
    t1 = datetime('now');
    dt = t1 - t0;
    fprintf('Startup time: %s',dt)

    fprintf('\n')
end

