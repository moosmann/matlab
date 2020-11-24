function nexuslog_name = pp_get_nexuslog_names( par )
% Helper function returning a cell of all paths to the nexus log files

% Scan log
nexuslog_name = {get_nexuslog_name( par.scan_path )};

% Addtional ref log
if ~isempty( par.ref_path )
    for nn = 1:numel( par.ref_path )
        rp = par.ref_path{nn};
        nln = get_nexuslog_name( rp );
        nexuslog_name = cat( 2, nexuslog_name, nln );
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nexuslog_name = get_nexuslog_name( logpath )

nexuslog_name = dir( sprintf('%s*_nexus.h5', logpath ) );
if numel( nexuslog_name ) == 1
    nexuslog_name = [nexuslog_name.folder filesep nexuslog_name.name];
end

