function [poolobj, poolsize] = OpenParpool( poolsize, use_cluster, tmp_folder, enforce_poolsize, poolsize_max)
% Open pool of parallel worker if it doesn't exist or if it is smaller than
% the desired pool size.
%
% poolsize: scalar. Number of workers desired.
% use_cluster : bool. use parpool on cluster with 256 cores using SLURM
%   scheduling.
% tmp_folder : string, directory with write permission, required for the
%   parpool to start
%
% Written by Julian Moosmann. Last modification: 2016-10-06, 2018-01-15
%
% [poolobj, poolsize] = OpenParpool(poolsize)

if nargin < 2
    use_cluster = 0;
end
if nargin < 3
    tmp_folder = '';
end
if nargin < 4
    enforce_poolsize = 0;
end
if nargin < 5
    poolsize_max = [];
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = toc;
cluster_poolsize = 250; % max 256

if isempty( tmp_folder )  
    tmp_folder = sprintf( '/gpfs/petra3/scratch/%s/', getenv( 'USER' ) );
end
CheckAndMakePath( tmp_folder )

% Desired poolsize
if poolsize < 1 && poolsize > 0
    numCores = feature('numCores');
    poolsize = max( floor( poolsize * numCores ), 1 );
end

% Max poolsize
if ~isempty( poolsize_max )
    poolsize = min( poolsize, poolsize_max );
end

% check if more than 1 worker is desired
if poolsize > 1
    
    % check if cluster computation is available
    hostname = getenv('HOSTNAME');    
    if use_cluster && sum( strcmp( hostname(1:8), {'max-disp','max-nova', 'max-wgs', 'max-wga'}) )
        % above nodes support SLURM scheduling
        clust = 'maxwell';
        poolsize = cluster_poolsize;
    else
        clust = 'local';
    end
    
    cpath = pwd;
    cd( tmp_folder );
    % check if already parpool exists
    if isempty( gcp('nocreate') )
        %if not open parpool
        fprintf('\n')
        poolobj = parpool( clust, poolsize);
    else
        % get current pool
        poolobj = gcp;
        
        % check if current pool is desired pool
        pflag = 0;
        if ~strcmp( poolobj.Cluster.Profile, clust )
            pflag = 1;
        end
        
        % check if number of workers of current pool is smaller than poolsize
        if strcmp( clust, 'local' )  && poolobj.NumWorkers < poolsize
            pflag = 1;
        end
        
        % if number of worker of current pool is larger than poolsize_max
        if strcmp( clust, 'local' )  && poolobj.NumWorkers > poolsize_max
            pflag = 1;
        end
        
        if enforce_poolsize
            pflag = 1;
        end
        
        % delete current pool and open new one
        if pflag
            fprintf('\n')
            poolobj.delete;
            poolobj = parpool( poolsize);
        end
    end
    
    cd( cpath );
end

if poolsize > 1
    fprintf( '\nParpool opened on %s using %u of %u workers', poolobj.Cluster.Profile, poolobj.NumWorkers, poolobj.Cluster.NumWorkers )
    fprintf( ' in %.1f s (%.2f min)', toc - t, ( toc - t ) / 60 )
end
