% CUDA_ART.m

clear all;

gpuDevice(3);

volSize = [1200 1200 1];
projSize = [1008 1];
phiRange = (0:180/200:(5795)/200*180)+50.4;

%per = randperm(length(phiRange));
%perTime = per+5400;
%per = mod(floor(((per-1)*(length(per)+1)/2)),length(per))+1;

%in_path = '/mnt/tomoraid3/user/haenschke/jakob/Odondomachus';
%in_prefix = 'ffc';

in_path = '/mnt/tomoraid3/user/haenschke/simulations/ARTcuda/time_series';
in_prefix = 'sino_1500fps.edf';
out_prefix = 'reco_tryOut_blocks3_1500fps';

projection_step_size = 50;
projs_per180deg = 200;

recOrig = volSize/2+[0 0 0];
origProj = [517 0.5 0.5];%projSize/2;

oversampling = 4;
osPos = -0.5+1/(2*oversampling)+(0:1/oversampling:(oversampling-1)/oversampling);

theta = 0;

sirt_interval = 100; % 50 next time!
iterations = 50; %50 next time!
start_iterations = 50;

gridSize = 32;
noThreads = 32;

k_proj = parallel.gpu.CUDAKernel('CUDA_ART_getProj_slice_2.ptx','CUDA_ART_getProj_slice_2.cu');
k_proj.ThreadBlockSize = [noThreads 1];
k_proj.GridSize = [gridSize 1];

k_corr = parallel.gpu.CUDAKernel('CUDA_ART_getCorr_slice_2.ptx','CUDA_ART_getCorr_slice_2.cu');
k_corr.ThreadBlockSize = [noThreads 1 1];
k_corr.GridSize = [gridSize 1];

R_th = [cosd(theta) 0 sind(theta); 0 1 0; -sind(theta) 0 cosd(theta)];

vol_g = parallel.gpu.GPUArray.zeros(volSize(1),volSize(2),volSize(3),'single');
vol_corr_g = parallel.gpu.GPUArray.zeros(volSize(1),volSize(2),volSize(3),'single');
%projs_temp = zeros(projSize(1),projSize(2),'single');

%y_slice = 1200;

%projs = zeros(projSize(1),length(phiRange),'single');
 
%projs_save = zeros(projSize(1),length(phiRange),'single');  --->!

% % load projections
% count_projs = 0;
% for p=1:length(phiRange)
%    count_projs = count_projs + 1;
%    inName = sprintf('%s/%s_%4.4i.edf', in_path, in_prefix, p);
%    projs_temp(:,:) = single(edfread(inName));
%    projs_save(:,count_projs) =  projs_temp(:,y_slice);
%    projs_temp(:,y_slice) = (projs_temp(:,y_slice)>0).*(projs_temp(:,y_slice)).*(projs_temp(:,y_slice)<=1) + (projs_temp(:,y_slice)>1)*1.0;
%    projs(:,count_projs) = -log(projs_temp(:,y_slice));
%    fprintf('projection no %d of %d loaded.\n',count_projs, length(phiRange));
% end

inName = sprintf('%s/%s', in_path, in_prefix);
projs = single(edfread(inName));
proj_g = gpuArray(squeeze(projs(:,:)));
proj_g_corr = parallel.gpu.GPUArray.zeros(projSize(1),length(phiRange),'single');
% imtool(projs);


 % start iterations
count = 0;
for it=1:start_iterations
tic;
start_i = 5795-projs_per180deg;
per_i = randperm(projs_per180deg);
for i=1:projs_per180deg
    pos = per_i(i)+start_i;
    count = count + 1;
    % proj_g_corr = 0*proj_g_corr;
    phi = phiRange(pos);    
    R = [cosd(phi) -sind(phi) 0; sind(phi) cosd(phi) 0; 0 0 1]*R_th;
    proj_g_corr = 0*proj_g_corr;
    
    for os=osPos
    % calculate the projection difference
    [proj_g_corr] = feval(k_proj,R,volSize,projSize,vol_g,recOrig,origProj,proj_g,pos,proj_g_corr,[os 0]);
    end
    
    proj_g_corr = proj_g_corr/(oversampling^2);
    
    % calculate the correction volume
    for os=osPos
    [vol_corr_g] = feval(k_corr,R,volSize,projSize,vol_corr_g,recOrig,origProj,proj_g_corr,pos,[0 0],[os 0]);
    [vol_corr_g] = feval(k_corr,R,volSize,projSize,vol_corr_g,recOrig,origProj,proj_g_corr,pos,[0 1],[os 0]);
    [vol_corr_g] = feval(k_corr,R,volSize,projSize,vol_corr_g,recOrig,origProj,proj_g_corr,pos,[1 0],[os 0]);
    [vol_corr_g] = feval(k_corr,R,volSize,projSize,vol_corr_g,recOrig,origProj,proj_g_corr,pos,[1 1],[os 0]);
    %fprintf('proj %d of %d (iteration %d of %d) done..\n',i,length(phiRange),it,iterations);
    end
    
    % calculate the correction
    if (mod(count,sirt_interval)==0)
       count = 0;
       vol_g = vol_g + vol_corr_g/sirt_interval;
       %vol_g = (vol_g>0).*vol_g;
       vol_corr_g = 0*vol_corr_g;
    end    
    %values(count) = proj_test(6,6);
end
if (count~=0)
    count = 0;
    vol_g = vol_g + vol_corr_g*count/sirt_interval;
    %vol_g = (vol_g>0).*vol_g;
    vol_corr_g = 0*vol_corr_g;
end
%vol_g = (vol_g>0).*vol_g;
fprintf('iteration %d of %d done.\n',it,start_iterations);
toc;
end

vol = gather(vol_g);
outName = sprintf('%s/%s_start.edf',in_path, out_prefix);
edfwrite(outName, double(vol), 'float32');
fprintf('start frame saved.\n');

count = 0;
save_count = 1;
for startProj = 5795-projs_per180deg-projection_step_size:-projection_step_size:1
for it=1:iterations
tic;
per_i = randperm(projs_per180deg);
perm_i = per_i - 1 + startProj;
for i=1:projs_per180deg
    count = count + 1;
    % proj_g_corr = 0*proj_g_corr;
    phi = phiRange(perm_i(i));    
    R = [cosd(phi) -sind(phi) 0; sind(phi) cosd(phi) 0; 0 0 1]*R_th;
    proj_g_corr = 0*proj_g_corr;
    
    for os=osPos
    % calculate the projection difference
    [proj_g_corr] = feval(k_proj,R,volSize,projSize,vol_g,recOrig,origProj,proj_g,perm_i(i),proj_g_corr,[os 0]);
    end
    
    proj_g_corr = proj_g_corr/(oversampling^2);
    
    % calculate the correction volume
    for os=osPos
    [vol_corr_g] = feval(k_corr,R,volSize,projSize,vol_corr_g,recOrig,origProj,proj_g_corr,perm_i(i),[0 0],[os 0]);
    [vol_corr_g] = feval(k_corr,R,volSize,projSize,vol_corr_g,recOrig,origProj,proj_g_corr,perm_i(i),[0 1],[os 0]);
    [vol_corr_g] = feval(k_corr,R,volSize,projSize,vol_corr_g,recOrig,origProj,proj_g_corr,perm_i(i),[1 0],[os 0]);
    [vol_corr_g] = feval(k_corr,R,volSize,projSize,vol_corr_g,recOrig,origProj,proj_g_corr,perm_i(i),[1 1],[os 0]);
    %fprintf('proj %d of %d (iteration %d of %d) done..\n',i,length(phiRange),it,iterations);
    end
    
    % calculate the correction
    if (mod(count,sirt_interval)==0)
       count = 0;
       vol_g = vol_g + vol_corr_g/sirt_interval;
       %vol_g = (vol_g>0).*vol_g;
       vol_corr_g = 0*vol_corr_g;
    end    
    %values(count) = proj_test(6,6);
end
if (count~=0)
    count = 0;
    vol_g = vol_g + vol_corr_g*count/sirt_interval;
    %vol_g = (vol_g>0).*vol_g;
    vol_corr_g = 0*vol_corr_g;
end
%vol_g = (vol_g>0).*vol_g;
fprintf('iteration %d of %d of block %d done.\n',it,iterations);
toc;
end
vol = gather(vol_g);
save_count = save_count + 1;
outName = sprintf('%s/%s_it%4.4i.edf',in_path, out_prefix, save_count-1);
edfwrite(outName, double(vol), 'float32');
fprintf('frame at block %d saved.\n',save_count-1);
end

%vol = gather(vol_g);
%outName = sprintf('%s/%s.edf',in_path, out_prefix);
%edfwrite(outName, double(vol), 'float32');
