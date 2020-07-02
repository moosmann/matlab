% CUDA_ART.m

clear all;

gpuDevice(2);

volSize = [100 100 100];
projSize = [752 672];
projSizeOrg = [752 672];
phiRange = 0:180/1200:(180-1/1200);

per = 1:length(phiRange);
per = mod(floor(((per-1)*(length(per)+1)/2)),length(per))+1;

in_path = '/home/moosmann/tomo/wildtype_05min_deadtime_05tomo_stage16p0_upwards_620mm_050ms_30p0keV/phase/TIElo_alpha4p02';
in_prefix = 'phase_'; 
out_path = '/home/moosmann/tomo/wildtype_05min_deadtime_05tomo_stage16p0_upwards_620mm_050ms_30p0keV/vol';

recOrig = volSize/2+[0 0 0];
origProj = [367.5 projSize(2)/2];

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
   inName = sprintf('%s/%s%04u.edf', in_path, in_prefix, p);
   im = single(pmedfread(inName));
   %im = im(:,floor(projSizeOrg(1)/2));
   projs(:,:,p) = im;
end

count = 0;
for it=1:iterations
tic;
fprintf('ITERATION %d OF %d.PROCESSING PROJECTION NUMBER OF %u PROJECTIONS:\n',it,iterations,length(phiRange));
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
    fprintf(' %4u',i);
    if mod(i,30)==0
        fprintf('\n')
    end
    
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

outName = sprintf('%s/rec_%s_%dx%dx%d.vol',out_path, in_prefix, volSize(1), volSize(2), volSize(3));
mexVolWrite(outName, double(vol), 'float32');