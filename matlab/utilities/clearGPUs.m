%% Free gpu memory occupied by parpool workers
g = gcp;
parfor n = 1:g.NumWorkers
    gpuDevice([]);
end