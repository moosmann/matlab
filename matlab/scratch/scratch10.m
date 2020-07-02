warning( 'off', 'parallel:gpu:device:DeviceDeprecated' );

gpu = parallel.gpu.GPUDeviceManager.instance;
tic;
num_gpu = gpuDeviceCount;
toc
for n = 1:num_gpu
    tic;
    m =  p.isAvailable(n);
    
    fprintf( '\n %u num %u, is available',n, m )
    t = toc;
end


%% parpool creation
fprintf( '\n Parpool creation and deleting times\n' )
tic
parpool( 'threads' );
t = toc;
fprintf( '\n parpool threads create : %f s\n', t )

pp = gcp( 'nocreate' );
tic
delete( pp );
t = toc;
fprintf( '\n parpool threads delete : %f s\n', t )

tic
parpool( 'local' );
t = toc;
fprintf( '\n parpool local create : %f s\n', t )

pp = gcp( 'nocreate' );
tic
delete( pp );
t = toc;
fprintf( '\n parpool local delete : %f s\n', t )


%% parpool local
%pp = gcp( 'nocreate' );

tic;
spmd(18,32)
    n = labindex;    
    %ngpu = gpuDevice;
    gpuDev  = gpuDevice;
    ngpu = gpuDev.Index;
    fprintf( ' labindex %u, gpu %u', n, ngpu )
end
toc

%% parpool threads
pp = gcp( 'nocreate' );
num_workers = pp.NumWorkers;
fprintf( '\n' )
parfor n = 1:num_workers
    ngpu = gpuDevice;
    fprintf( ' worker %u, gpu %u', n, ngpu )
end