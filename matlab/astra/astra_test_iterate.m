tic
aclear;
angleInc = 1;
angularIncrementInRad = pi/500;
maxIter = 500;
initializer = 0;
recType = 'SIRT_CUDA';
%recType = 'CGLS_CUDA';

% Read sinogramm
ParentPath = '/export/scratch1/moosmann/APS_32ID-C_InVivoImaging_GUP34879_2013-07-30/vol/Xenopus_inVivo/Jul29_15-10_urea_stage27p0_30p0keV_0700mm_15ms_0500proj_scantime20s_deadtime8min/tomo01';
%sino = Readstack(sprintf('%s/sino_intProjInc1__0499x1966x0100',ParentPath));
sino = imread(sprintf('%s/sino_intProjInc1__0499x1966x0100/slice_0050.tif',ParentPath));
% subtract mean
sino = sino -mean(sino(:));

% parameter
[NumProjFull, dimHor] = size(sino);
projUsed = 1:angleInc:NumProjFull;
NumProjRed = numel(projUsed);
dimHorVol = dimHor;
%dimHorVol = 2*round(dimHor/sqrt(2)/2);
angles = angularIncrementInRad*((1:NumProjFull)-1);

%% Geometry, id's, cfg
% Creat volume and projection geometry
vol_geom = astra_create_vol_geom([dimHorVol,dimHorVol]);
proj_geom = astra_create_proj_geom('parallel',1.0,size(sino,2),angles(projUsed));

% Create the sinogram data object
sino_id = astra_mex_data2d('create', '-sino', proj_geom, sino);

% Create a data object for the reconstruction

rec_id = astra_mex_data2d('create', '-vol', vol_geom,initializer);
cfg = astra_struct(recType);
% Parameters for a reconstruction algorithm using the GPU
cfg.ReconstructionDataId = rec_id;
cfg.ProjectionDataId = sino_id;
cfg.Options.GPUindex = 2; % Use GPU #1 for the reconstruction. (The default is #0.)
% Create the algorithm object from the configuration structure
alg_id = astra_mex_algorithm('create', cfg);

% Preallocation
rec = zeros([dimHorVol, dimHorVol, maxIter],'single');
res = zeros(maxIter,'single');

%% Iterate the algorithm one at a time, keeping track of errors
for nn = 1:maxIter
    PrintNum(nn)
    % Run a single iteration
    astra_mex_algorithm('iterate', alg_id, 1);
    % Get and save result
    rec(:,:,nn) = astra_mex_data2d('get_single', rec_id);
    %astra_mex_data2d('set',rec_id,FilterMedian(rec(:,:,nn)));

    res(nn) = astra_mex_algorithm('get_res_norm', alg_id);
    
end

% Clean up. Note that GPU memory is tied up in the algorithm object,
% and main RAM in the data objects.
astra_mex_algorithm('delete', alg_id);
astra_mex_data2d('delete', rec_id);
astra_mex_data2d('delete', sino_id);
