function OpenParpool(poolsize)
% Open pool of parallel worker if it doesn't exist or if it is smaller than
% the desired pool size.
%
% poolsize: scalar. Number of workers desired.


%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check if more than 1 worker is desired
if poolsize > 1
    
    % check if already parpool exists
    if isempty(gcp('nocreate'))
        %if not open parpool
        parpool( poolsize);
    else
        % get current pool
        poolobj = gcp;
        % if number of workers of current pool is smaller than poolsize,
        % delete current pool and open new one
        if poolobj.NumWorkers < poolsize
            fprintf('\n')
            poolobj.delete;
            parpool( poolsize);
        end        
    end
end
