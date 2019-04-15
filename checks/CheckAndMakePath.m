function CheckAndMakePath(name, deleteFiles, beamtimeID_regexp)
% Check if path 'name' exists and if not create it.

if nargin < 2
    deleteFiles = 0;
end
if nargin < 3
    beamtimeID_regexp = '';
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist(name,'dir')
    mkdir(name);
else
    if deleteFiles
        % Safety check
        
        expr = '11001978|c20160920_000_diana|c20160913_000_synload|11003950|11004016|11003288|11003440|11004263|11004936|11005553';
        if ~isempty( beamtimeID_regexp )
            expr = [expr '|' beamtimeID_regexp];
        end
        
        if isempty( expr )
            return
        end
        t = regexp( name, expr, 'once' );
        if ~isempty( t )
            s = dir( [name filesep '*.tif'] );
            parfor nn = 1:numel( s )
                filename = [s(nn).folder filesep s(nn).name ];
                if ~s(nn).isdir
                    delete( filename )
                end
            end
            
        end

    end
end
