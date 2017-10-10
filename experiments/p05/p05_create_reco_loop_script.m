function p05_create_reco_loop_script( raw_path, folder_pattern)
% Create template file to loop over reconstruction for all data sets, i.e.
% folder, found in 'raw_path' matching the 'folder_pattern'. Created file
% will be immediately opened in the MATLAB editor for editing.
%
% raw_path : string. Default: use present working directory. path to scan
%   for data sets. 
% folder_pattern : string. Default: ''. only add folders matching pattern,
%   e.g. 'dataSetNamePrefix*'. Asterisk (*) is required to match pattern.
%   
% Written by Julian Moosmann, 2017-10-10. Last version: 2017-10-10
%
% p05_create_reco_loop_script( raw_path, folder_pattern)

% TODO: Issue: where to create file
% TODO: Issue: related to file location: add search path

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1   
    raw_path = pwd;
end
if nargin < 2
    folder_pattern = '';
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% directories
[raw_path_prefix, name]= fileparts( raw_path );
if ~strcmp( name, 'raw' )
    fprintf( '\n WARNING: Name of current directory (%s) is not ''raw''', name)
end
proc_path = [raw_path_prefix  filesep 'processed'];

%% Folder
folders = dir( [raw_path filesep folder_pattern]);
isub = [folders(:).isdir];
folders = {folders(isub).name};
fprintf( '\n Found %u data sets in ''%s'' matching the pattern:', numel(folders), raw_path )
for nn = 1:numel( folders )
    fprintf( '\n %4u : %s', nn, folders{nn})
end

%% Create new, non-existing filename
nn = 0;
filename_prefix = [proc_path filesep 'p05_reco_loop_script'];
filename = sprintf( '%s_%02u.m', filename_prefix, nn );
while exist( filename, 'file')
    nn = nn + 1;
    filename = sprintf( '%s_%02u.m', filename_prefix, nn );
end

%% Create template %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read top and bottom parts from files
p = '/asap3/petra3/gpfs/common/p05/jm/matlab/experiments/p05/';
top = fileread( sprintf('%s%s', p, 'reco_loop_template_top.m' ) );
bottom = fileread( sprintf( '%s%s', p, 'reco_loop_template_bottom.m' ) );

% Write top: function defition, help arguments, etc
fid = fopen( filename, 'w');
func_name = sprintf( 'p05_reco_loop_script_%02u', nn);
fprintf(fid, 'function %s%s', func_name, top);

% Write Data Sets
fprintf(fid, '\nraw_path = ''%s/'';', raw_path);
for nn = 1:numel( folders )
    fprintf(fid, '\n' );
    fprintf(fid, '\nscan_path = [raw_path %s]; ADD', folders{nn});    
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
