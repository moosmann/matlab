function p05_create_reco_loop_script( raw_path, folder_pattern, appendix)
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
% appendix : str. Default: ''. String to append to filename of the script.
%   If empty a two-digit number is used to avoid overwriting existing
%   scripts.
%   
% Written by Julian Moosmann, 2017-10-10. Last version: 2017-10-14
%
% p05_create_reco_loop_script( raw_path, folder_pattern, appendix)

% TODO: How to handle parameters!!!!
% TODO: Option to overwrite existing file 
% TODO: Where to create file
% TODO: Related to file location: add search path
% TODO: Add creation time and user name to template

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
folders(ismember(folders,{'.','..'})) = [];
for nn = 1:numel( folders )
    fprintf( '\n %4u : %s', nn, folders{nn})
end

%% Create new, non-existing filename

func_name = 'p05_reco_loop_script';
parent_path = [proc_path filesep];
if isempty( appendix )
    nn = 0;
    func_name = sprintf( '%s_%02u', func_name, nn);
    filename = sprintf( '%s_%s.m', parent_path, func_name );
    while exist( filename, 'file')
        nn = nn + 1;
        func_name = sprintf( '%s_%02u', func_name, nn);
        filename = sprintf( '%s_%s.m', parent_path, func_name );        
    end
    
else
    func_name = sprintf( '%s_%s', func_name, appendix);
    filename = sprintf( '%s%s.m', parent_path, func_name);
end


%% Create template %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read top and bottom parts from files
p = '/asap3/petra3/gpfs/common/p05/jm/matlab/experiments/p05/';
top = fileread( sprintf('%s%s', p, 'reco_loop_template_top.m' ) );
bottom = fileread( sprintf( '%s%s', p, 'reco_loop_template_bottom.m' ) );

% Write top: function defition, help arguments, etc
fid = fopen( filename, 'w');
func_name = sprintf( '%s', func_name);
fprintf(fid, 'function %s%s', func_name, top);

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
