% Test parpool creation and combination of parpool and spmd with GPUs.
warning( 'off', 'parallel:gpu:device:DeviceDeprecated' );

gpu = parallel.gpu.GPUDeviceManager.instance;
tic;
num_gpu = gpuDeviceCount;
toc
fprintf( 'GPU device count: %u', num_gpu)
for n = 1:num_gpu
    tic;
    m =  gpu.isAvailable(n);
    
    fprintf( '\n %u num %u, is available',n, m )
    t = toc;
end


%% parpool creation
fprintf( '\n Parpool creation and deleting times\n' )

%% parpool threads
tic
parpool( 'threads' );
t = toc;
fprintf( '\n parpool threads create : %f s\n', t )

pp = gcp( 'nocreate' );
tic
delete( pp );
t = toc;
fprintf( '\n parpool threads delete : %f s\n', t )

%% parpool local
tic
parpool( 'local' );
t = toc;
fprintf( '\n parpool local create : %f s\n', t )

pp = gcp( 'nocreate' );
tic
delete( pp );
t = toc;
fprintf( '\n parpool local delete : %f s\n', t )


%% Parpool and GPUs
pp = gcp( 'nocreate' );
num_workers = pp.NumWorkers;
fprintf( '\n' )
fprintf('\nparpool and GPU device' )
parfor n = 1:num_workers
    gpuDev = gpuDevice;
    ngpu = gpuDev.Index;
    fprintf( '\n  worker %u, gpu %u', n, ngpu )
end

%% SPMD

fprintf( '\nSPMD')
%pp = gcp( 'nocreate' );
tic;
N = pp.NumWorkers;
fprintf('\n numWorkers: %u', N );
fprintf('\n numGPU> %u', num_gpu);
fprintf('\n')
spmd(num_gpu, N)
    n = labindex;        
    gpuDev  = gpuDevice;
    ngpu = gpuDev.Index;
    fprintf( ' labindex %u, gpu %u', n, ngpu )
end
toc