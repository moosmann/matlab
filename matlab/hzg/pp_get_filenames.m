function [proj_names, ref_names, ref_full_path, dark_names, par] =  pp_get_filenames( par )
% Get filenames from disk

scan_path = par.scan_path;
ref_path = par.ref_path;

% Projection file names
proj_names = FilenameCell( [scan_path, '*.img'] );
par.raw_data = 0;
if isempty( proj_names )
    proj_names =  FilenameCell( [scan_path, '*img*.tif'] );
    par.raw_data = 0;
end
if isempty( proj_names )
    proj_names =  FilenameCell( [scan_path, '*img*.raw'] );
    par.raw_data = 1;
end
if isempty( proj_names )
    proj_names =  FilenameCell( [scan_path, '*proj*.tif'] );
    par.raw_data = 0;
end
if isempty( proj_names )
    proj_names =  FilenameCell( [scan_path, '*proj*.raw'] );
    par.raw_data = 1;
end
% Ref file names
[ref_names, ref_full_path] = get_ref_names( scan_path );

% Dark file names
dark_names = FilenameCell( [scan_path, '*.dar'] );
if isempty( dark_names )
    dark_names = FilenameCell( [scan_path, '*dar.tif'] );
end
if isempty( dark_names )
    dark_names =  FilenameCell( [scan_path, '*dar*.raw'] );
end
if isempty( dark_names )
    dark_names = FilenameCell( [scan_path, '*dar*.tif'] );
end

%% Add additional ref_paths
if ~isempty( ref_path )
    % Backup old parameter
    par.ref_path_old = par.ref_path;
    % In case non-cell parameter is given
    if ~iscell( ref_path)
        ref_path = {ref_path};
    end
    % Preallocation
    ref_path_delete = zeros( size( ref_path ), 'logical' );
    ref_path_to_add = {};
    % Loop over ref_path entries
    for n = 1:numel( ref_path )
        rp = ref_path{n};
        % Check for relative/absolute path
        if rp(1) ~= filesep
            % Full path
            rp_full = [par.raw_path filesep rp];
        else
            rp_full = rp;
        end
        
        % Check for asterisk indirectly
        if isfolder( rp_full )
            CheckTrailingSlash(rp_full)
            ref_path{n} = rp_full;
        else
            % To be removed entry
            ref_path_delete(n) = 1;
            % Search for matching folders
            RemoveTrailingSlash( rp_full )
            d = dir( rp_full );
            if ~isempty( dir )
                % Loop over matching patterns
                for m = 1:numel( d )
                    % Full path
                    rp_full = [d(m).folder filesep d(m).name filesep];
                    % Add path to cell of ref paths
                    if isfolder( rp_full ) && ~strcmpi( rp_full, scan_path ) && ~sum( contains( ref_path_to_add, rp_full) )
                        ref_path_to_add = cat( 2, ref_path_to_add, rp_full );
                    end
                end
            end
        end
        
    end
    
    % Remove entries
    ref_path( ref_path_delete ) = [];
    ref_path = cat( 2, ref_path, ref_path_to_add );
    par.ref_path = ref_path;
    fprintf( '\n Additional reference paths included:')
    fprintf( '\n %s', ref_path{:})
end

if ~isempty( ref_path )
    % Loop over ref paths
    for n = 1:numel( ref_path )
        [rn, rfp] = get_ref_names( ref_path{n} );
        ref_names = cat( 2, ref_names, rn );
        ref_full_path = cat( 2, ref_full_path, rfp );
    end
    
end

%% Hack due to rewriting of tomoscan_flikit
if isempty( ref_names )
    h5log = dir( sprintf('%s*_nexus.h5', scan_path) );
    h5log = [h5log.folder filesep h5log.name];
    % images
    stimg_name.value = unique( h5read( h5log, '/entry/scan/data/image_file/value') );
    stimg_name.time = h5read( h5log,'/entry/scan/data/image_file/time');
    stimg_key.value = h5read( h5log,'/entry/scan/data/image_key/value');
    stimg_key.time = double( h5read( h5log,'/entry/scan/data/image_key/time') );
    
    % File names
    proj_names = stimg_name.value(stimg_key.value==0)';
    ref_names = stimg_name.value(stimg_key.value==1)';
    dark_names = stimg_name.value(stimg_key.value == 2 )';
    
    ref_full_path = cellfun( @(a) [scan_path filesep a], ref_names, 'UniformOutput', 0 );
    
    %% name check for petra iii current is now useless    
    cprintf( 'Red', '\nWarning: HACK. USELESS CROSS CHECK OF FILENAMES AND LOGFILE' )
end

function [ref_names, ref_full_path] = get_ref_names( scan_path )
[ref_names, ref_full_path] = FilenameCell( [scan_path, '*.ref'] );
if isempty( ref_names )
    [ref_names, ref_full_path]  = FilenameCell( [scan_path, '*ref.tif'] );
end
if isempty( ref_names )
    [ref_names, ref_full_path]  =  FilenameCell( [scan_path, '*flat*.tif'] );
end
if isempty( ref_names )
    [ref_names, ref_full_path]  =  FilenameCell( [scan_path, '*ref*.raw'] );
end
if isempty( ref_names )
    [ref_names, ref_full_path]  =  FilenameCell( [scan_path, '*flat*.raw'] );
end
if isempty( ref_names )
    [ref_names, ref_full_path]  = FilenameCell( [scan_path, '*ref*.tif'] );
end