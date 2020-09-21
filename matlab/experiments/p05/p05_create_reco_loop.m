function p05_create_reco_loop( raw_path, scan_name_pattern, out_path, name)
% Create template script to loop reconstruction over all data sets i.e.
% folders found under 'raw_path' matching the 'scan_name_pattern'. The
% created script will be opened immediately in the MATLAB editor for
% further editing. 
%
% !! Uses parameter setting of the main reco script 'p05_reco' if not
% changed in the created templated file. See help within created template!!
%
% ARGUMTENTS:
% raw_path : string. Default or if empty: use present working directory, path to scan
%   for data sets. 
% scan_name_pattern : string. Default: ''. only add folders matching pattern,
%   e.g. 'dataSetNamePrefix*'. Asterisk (*) is required to match pattern.
% out_path : path where loop script will be saved. It's recommended to a
%   use path within the MATLAB search path. Default:
%   $MATLAB_SEARCH_PATH/experiments/p05/data/$USER/
% name : str. Default: ''. String to append to filename of the script.
%   A 3-digit running index is used to avoid overwriting existing scripts.
%   
% Written by Julian Moosmann, 2017-10-10. Last version: 2018-07-22
%
% p05_create_reco_loop( raw_path, scan_name_pattern, name, out_path)

% TODO: Option to overwrite existing file 
% TODO: file location: add to search path

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1   
    raw_path = pwd;
end
if nargin < 2
    scan_name_pattern = '';
end
if nargin < 3
    out_path = '';
end
if nargin < 4
    name = '';
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum( strcmp( raw_path, {'.',''}) )
    raw_path = pwd;
end

%% directories
[pathstr, raw_name]= fileparts( raw_path );
[~, beamtimeid] = fileparts( pathstr );
if isnan( str2double( beamtimeid ) )
    beamtimeid = '';
end
if ~strcmp( raw_name, 'raw' )
    fprintf( '\n WARNING: Name of current directory (%s) is not ''raw''', raw_name)
    input( '\n If path is correct type ENTER to continue, otherwise abort.\n' );    
end
if isempty( out_path )    
    out_path = [userpath '/experiments/p05/data/' getenv('USER') filesep];
end
CheckTrailingSlash( out_path )
CheckAndMakePath( out_path )

%% Folders to read
folders = dir( [raw_path filesep scan_name_pattern]);
isub = [folders(:).isdir];
folders = {folders(isub).name};
fprintf( '\n Found %u data sets in ''%s'' matching the pattern:', numel(folders), raw_path )
folders(ismember(folders,{'.','..'})) = [];
for nn = 1:numel( folders )
    fprintf( '\n %4u : %s', nn, folders{nn})
end

%% Create new, non-existing filename

func_name0 = 'p05_reco_loop_';
parent_path = out_path;
if ~isempty( name )
    func_name0 = [func_name0 name '_' ];
end
func_name0 = [ func_name0 beamtimeid ];
nn = 0;
func_name = sprintf( '%s_%03u', func_name0, nn);
filename = sprintf( '%s%s.m', parent_path, func_name );
while exist( filename, 'file')
    nn = nn + 1;
    func_name = sprintf( '%s_%03u', func_name0, nn);
    filename = sprintf( '%s%s.m', parent_path, func_name );
end

%% Create template %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read top and bottom parts from files
p = [userpath '/matlab/experiments/p05/'];
top = fileread( sprintf('%s%s', p, 'reco_loop_template_top.m' ) );
middle = fileread( sprintf('%s%s', p, 'reco_loop_template_middle.m' ) );
middle_2 = fileread( sprintf('%s%s', p, 'reco_loop_template_middle_2.m' ) );
bottom = fileread( sprintf( '%s%s', p, 'reco_loop_template_bottom.m' ) );

% Write top: function defition, help arguments
fid = fopen( filename, 'w');
func_name = sprintf( '%s', func_name);
fprintf(fid, 'function %s%s', func_name, top);
fprintf(fid, 'Created on %s by %s\n', date, getenv('USER') );

% Write middle: default arguments
fprintf(fid, '\n%s', middle);

% Paramter set
fidrec = fopen( [userpath '/matlab/experiments/p05/p05_reco.m'], 'r' );
writetag = 0;
while 1
    c = textscan(fidrec,'%s',1, 'Delimiter', {'\n'});
    c = c{1}{1};
    if strcmp( c, '%% PARAMETERS / SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        writetag = 1;
        continue
    end    
    if strcmp( c, '%%% END OF PARAMETERS / SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')        
        break        
    end
    if writetag
        fprintf( fid, '\n%s', c );
    end   
end
fclose( fidrec );

% Write 2nd middle parte
fprintf(fid, '\n%s', middle_2);

% Write Data Sets
fprintf(fid, '\nraw_path = ''%s/'';', raw_path);
for nn = 1:numel( folders )
    fprintf(fid, '\n' );
    fprintf(fid, '\nscan_path = [raw_path ''%s'']; ADD', folders{nn});    
end

% Write bottom: call reco loop
fprintf(fid, '\n%s', bottom);
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\n Reco loop template script created:\n %s', filename )
fprintf( '\n Opening template script in MATLAB editor.' )
fprintf( '\n Follow instructions in template script to start looping.' )
edit( filename )
fprintf( '\n' )
