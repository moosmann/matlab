% CUDA_ART.m

clear all;

gpuDevice(2);

volSize = [141 141 141];
projSize = [201 201];
phiRange = 0:1:179;

per = 1:length(phiRange);
per = mod(floor(((per-1)*(length(per)+1)/2)),length(per))+1;

in_path = '/mnt/tomoraid3/user/haenschke/simulations/ARTcuda/test_ART';
in_prefix = 'test_ART'; 

recOrig = volSize/2+[0 0 0];
origProj = projSize/2;

theta = 0;

sirt_interval = 90;
iterations = 20;

k_proj = parallel.gpu.CUDAKernel('CUDA_ART_getProj.ptx','CUDA_ART_getProj.cu');
k_proj.ThreadBlockSize = [64 1 1];
k_proj.GridSize = [64 1];

k_corr = parallel.gpu.CUDAKernel('CUDA_ART_getCorr.ptx','CUDA_ART_getCorr.cu');
k_corr.ThreadBlockSize = [64 1 1];
k_corr.GridSize = [64 1];

R_th = [cosd(theta) 0 sind(theta); 0 1 0; -sind(theta) 0 cosd(theta)];

vol_g = parallel.gpu.GPUArray.zeros(volSize(1),volSize(2),volSize(3),'single');
vol_corr_g = parallel.gpu.GPUArray.zeros(volSize(1),volSize(2),volSize(3),'single');
projs = zeros(projSize(1),projSize(2),length(phiRange),'single');

% load projections
for p=1:length(phiRange)
   inName = sprintf('%s/%s_%4.4i.edf', in_path, in_prefix, p-1);
   projs(:,:,p) = single(edfread(inName));
end

count = 0;
for it=1:iterations
tic;
for i=1:length(phiRange)
    count = count + 1;
    
    proj_g = gpuArray(squeeze(projs(:,:,i)));
    phi = phiRange(per(i));    
    R = [cosd(phi) -sind(phi) 0; sind(phi) cosd(phi) 0; 0 0 1]*R_th;
    
    % calculate the projection difference
    [proj_g] = feval(k_proj,R,volSize,projSize,vol_g,recOrig,origProj,proj_g);
    
    % calculate the correction volume    
    [vol_corr_g] = feval(k_corr,R,volSize,projSize,vol_corr_g,recOrig,origProj,proj_g,[0 0]);
    [vol_corr_g] = feval(k_corr,R,volSize,projSize,vol_corr_g,recOrig,origProj,proj_g,[0 1]);
    [vol_corr_g] = feval(k_corr,R,volSize,projSize,vol_corr_g,recOrig,origProj,proj_g,[1 0]);
    [vol_corr_g] = feval(k_corr,R,volSize,projSize,vol_corr_g,recOrig,origProj,proj_g,[1 1]);
    fprintf('proj %d of %d (iteration %d of %d) done..\n',i,length(phiRange),it,iterations);
    
    % calculate the correction
    if (mod(count,sirt_interval)==0)
       count = 0;
       vol_g = vol_g + vol_corr_g/sirt_interval;
       vol_g = (vol_g>0).*vol_g;
       vol_corr_g = 0*vol_corr_g;
    end    
    %values(count) = proj_test(6,6);
end
if (count~=0)
    vol_g = vol_g + vol_corr_g*count/sirt_interval;
    vol_g = (vol_g>0).*vol_g;
    vol_corr_g = 0*vol_corr_g;
end
toc;
end

vol = gather(vol_g);

outName = sprintf('%s/rec_%s_%dx%dx%d.vol',in_path, in_prefix, volSize(1), volSize(2), volSize(3));
mexVolWrite(outName, double(vol), 'float32');