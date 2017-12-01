function p05_create_reco_loop_script( raw_path, folder_pattern, appendix, out_path)
% Create template file to loop over reconstruction for all data sets, i.e.
% folder, found in 'raw_path' matching the 'folder_pattern'. Created file
% will be immediately opened in the MATLAB editor for editing.
%
% !! Uses parameter setting of the main reco script 'p05_reco' if not
% changed in the created templated file. See help within created template!!
%
% ARGUMTENTS:
% raw_path : string. Default: use present working directory. path to scan
%   for data sets. 
% folder_pattern : string. Default: ''. only add folders matching pattern,
%   e.g. 'dataSetNamePrefix*'. Asterisk (*) is required to match pattern.
% out_path : path where loop script will be stored. It's recommended to a
%   path within the MATLAB search path.
% appendix : str. Default: ''. String to append to filename of the script.
%   If empty a two-digit number is used to avoid overwriting existing
%   scripts.
%   
% Written by Julian Moosmann, 2017-10-10. Last version: 2017-12-01
%
% p05_create_reco_loop_script( raw_path, folder_pattern, appendix, out_path)

% TODO: Option to overwrite existing file 
% TODO: file location: add search path

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1   
    raw_path = pwd;
end
if nargin < 2
    folder_pattern = '';
end
if nargin < 3
    appendix = '';
end
if nargin < 4
    out_path = '';
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp( raw_path, '.')
    raw_path = pwd;
end

%% directories
[~, name]= fileparts( raw_path );
%[raw_path_prefix, name]= fileparts( raw_path );
if ~strcmp( name, 'raw' )
    fprintf( '\n WARNING: Name of current directory (%s) is not ''raw''', name)
end
if isempty( out_path )    
    out_path = [userpath '/experiments/p05/data/' getenv('USER') filesep];
end
CheckTrailingSlash( out_path )
CheckAndMakePath( out_path )

%% Folders to read
folders = dir( [raw_path filesep folder_pattern]);
isub = [folders(:).isdir];
folders = {folders(isub).name};
fprintf( '\n Found %u data sets in ''%s'' matching the pattern:', numel(folders), raw_path )
folders(ismember(folders,{'.','..'})) = [];
for nn = 1:numel( folders )
    fprintf( '\n %4u : %s', nn, folders{nn})
end

%% Create new, non-existing filename

func_name0 = 'p05_reco_loop_script';
%processed_path = [raw_path_prefix  filesep 'processed'];
%parent_path = [processed_path filesep];
parent_path = out_path;

if ~isempty( appendix )
    func_name0 = [func_name0 '_', appendix ];
end
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
p = '/asap3/petra3/gpfs/common/p05/jm/matlab/experiments/p05/';
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
fidrec = fopen( '/asap3/petra3/gpfs/common/p05/jm/matlab/experiments/p05/p05_reco.m', 'r' );
writetag = 0;
while 1
    c = textscan(fidrec,'%s',1, 'Delimiter', {'\n'});
    c = c{1}{1};
    if strcmp( c, '%% PARAMETERS / SETTINGS %%')
        writetag = 1;
        continue
    end    
    if strcmp( c,  '%% END OF PARAMETERS / SETTINGS %%')        
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
