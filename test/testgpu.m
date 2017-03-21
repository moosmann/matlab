x=4097;
np=1;
sl=1;

mem = (x*x + x*np)*sl*4/1024^3;

fprintf( '\nmemory estimate for reco, single precision: %g GB', mem)

gpu = gpuDevice;
fprintf( '\ngpu memory: %g GP', gpu.AvailableMemory/1024^3)

v = astra_parallel3D(ones(x,np,sl,'single'),(0:np-1)/np*pi,0,[x x sl]);

fprintf( '\nvolume: %g GB', Bytes(v, 3))

fprintf( '\n')