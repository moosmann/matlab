function testAstra(varargin)
% Run FBP and ART reconstruction for a fiven set of parameters. 
%
% SART_CUDA only works with masking of the volume. EM_CUDA is not working.
%
% Written by Julian Moosmann, last version 2013-10-28

tic
aclear;
%% Parse input arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addParamValue(p,'Experiment','Xenopus4cell',@ischar);
addParamValue(p,'SubsetIndex',12,@isscalar);
addParamValue(p,'angleInc',1,@isscalar);
addParamValue(p,'horInc',1,@isscalar);
addParamValue(p,'filtLowFreq',0,@isscalar);
addParamValue(p,'maxIter',150,@isscalar);
addParamValue(p,'minIter',50,@isscalar);
addParamValue(p,'recType','SART_CUDA',@(x) ischar(x) || isscalar(x));
addParamValue(p,'initializer','',@(x) ischar(x) || isscalar(x));
addParamValue(p,'maskDiskDecr',8,@ischar);
addParamValue(p,'ParentPath','/export/scratch1/moosmann/art/art_vs_fbp/',@ischar);
addParamValue(p,'showFigs',1,@isscalar);
parse(p,varargin{:})
% Print parameters
hyphens(1:100) = '-';
fprintf('%s\nParamters:\n',hyphens)
disp(p.Results)
%% Set parameter
Experiment = p.Results.Experiment;
SubsetIndex = p.Results.SubsetIndex;
angularIncrementInRad = 2*pi/1599;
angleInc = p.Results.angleInc;
horInc = p.Results.horInc;
filtLowFreq = p.Results.filtLowFreq;
minIter = p.Results.minIter;
maxIter = p.Results.maxIter;
if ischar(p.Results.recType)
    recType = p.Results.recType;
elseif isscalar(p.Results.recType)    
    algorithm{1} = 'SIRT_CUDA';
    algorithm{2} = 'SART_CUDA';
    algorithm{3} = 'CGLS_CUDA';
    algorithm{4} = 'EM_CUDA';
    recType = algorithm{2};
end
initializer = p.Results.initializer;
maskDiskDecr = p.Results.maskDiskDecr;
ParentPath = p.Results.ParentPath;
showFigs = p.Results.showFigs;
%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read a sinogram
[sino, Experiment, sinoStr] = getSino(Experiment,SubsetIndex);
%sino = sino - repmat(mean(sino,2),[1 size(sino,2)]);
sino = sino -mean(sino(:));
fprintf('%15s: ''%s''\n\n','Subset',sinoStr)
% Make output folder
CheckTrailingSlash(ParentPath)
% Downsample number of pixels
sino = sino(:,1:horInc:end);
[NumProjFull, dimHor] = size(sino);
projUsed = 1:angleInc:NumProjFull;
NumProjRed = numel(projUsed);
dimHorVol = dimHor;
%dimHorVol = 2*round(dimHor/sqrt(2)/2);
% roi
if dimHorVol < dimHor
    x = 1:dimHorVol;
else
    dx = round(dimHorVol*(2-sqrt(2))/4*1.1);
    x = dx:dimHorVol-dx;
end
angles = angularIncrementInRad*((1:NumProj)-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Full FBP
fprintf(' Reconstruction:  FBP_CUDA, %4u of %u projections. ',NumProjFull,NumProjFull)
% Save sino
outputPath = MakePath('%s%s/%s/sino_pixInc%02u_projInc%02u_projUsed%04uof%04u',ParentPath,Experiment,sinoStr,horInc,angleInc,NumProjRed,NumProjFull);
filename = sprintf('%ssino_proj%04u_pix%04u',outputPath,size(sino));
WriteImage(filename,sino,{'tif','png'})
% Creat volume and projection geometry
vol_geom = astra_create_vol_geom([dimHorVol,dimHorVol]);
proj_geom = astra_create_proj_geom('parallel',1.0,size(sino,2),angles);
% Create the sinogram data object
sino_id = astra_mex_data2d('create', '-sino', proj_geom, sino);
% Create reconstruction
[fbp_id, fbp.full] = test_astra_fbp_cuda(proj_geom,vol_geom,sino_id);
astra_mex_data2d('delete', fbp_id);
astra_mex_data2d('delete', sino_id);
fprintf('Elapsed time: %g s\n',toc)

%% Reduced FBP
% Downsample number of projections
sino = sino(projUsed,:);
sino = sino -mean(sino(:));
fprintf(' Reconstruction:  FBP_CUDA, %4u of %u projections.',NumProjRed,NumProjFull)
% Creat volume and projection geometry
vol_geom = astra_create_vol_geom([dimHorVol,dimHorVol]);
proj_geom = astra_create_proj_geom('parallel',1.0,size(sino,2),angles(projUsed));
% Create the sinogram data object
sino_id = astra_mex_data2d('create', '-sino', proj_geom, sino);
% Create reconstruction
[fbp_id, fbp.red] = test_astra_fbp_cuda(proj_geom,vol_geom,sino_id);
astra_mex_data2d('delete', fbp_id);
fbp.blur = FilterBlur(fbp.red,[3 3],2);
% Save sino
filename = sprintf('%ssino_proj%04u_pix%04u',outputPath,size(sino));
WriteImage(filename,sino,{'tif','png'})
fprintf(' Elapsed time: %g s\n',toc)

%% ART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(' Reconstruction: %9s, %4u of %u projections. Iteration step:\n  ',recType,NumProjRed,NumProjFull)
% Filter low requencies
if filtLowFreq > 0
    [sino, filtLowFreqStr] = FilterLowFreq(sino,40,[1 120]);
    filename = sprintf('%s/sino_proj%04u_pix%04u_%s',outputPath,size(sino),filtLowFreqStr);
    WriteImage(filename,sino,{'tif','png'})
    %sino = FilterLowFreq(sino,9,[1 21]);
    % Save filtered sino
    outputPath = sprintf('%sreco',outputPath);
    outputPath = sprintf('%s_%s',outputPath,filtLowFreqStr);
else
    outputPath = sprintf('%sreco',outputPath);
end
% Initializer for ART
if ~isempty(initializer)
    switch lower(initializer)
        case 'fbp'
            outputPath = sprintf('%s_initializer%s',outputPath,initializer);
            initializer = fbp.full;
        case 'fbpblur'
            outputPath = sprintf('%s_initializer%s',outputPath,initializer);
            initializer = fbp.blur;
        case 'fbpcenblur'
            outputPath = sprintf('%s_initializer%s',outputPath,initializer);
            initializer = fbp.red;
            cen = [1:round(0.4*dimHorVol) round(0.6*dimHorVol):dimHorVol];
            initializer(cen,:) = 0;
            initializer(:,cen) = 0;
            initializer = FilterBlur(initializer);            
    end
            
else
    initializer = 0;
end
% Create a data object for the reconstruction
rec_id = astra_mex_data2d('create', '-vol', vol_geom,initializer);
cfg = astra_struct(recType);
if maskDiskDecr > 0
    mask = meshgrid((1:dimHorVol)-dimHorVol/2-0.5,(1:dimHorVol)-dimHorVol/2-0.5);
    mask = single(sqrt(mask.^2+mask'.^2) < dimHorVol/2 -maskDiskDecr);
    % Create a data object for the mask
    mask_id = astra_mex_data2d('create', '-vol', vol_geom, mask);
    cfg.option.ReconstructionMaskId = mask_id;
    outputPath = sprintf('%s_maskDecr%02u',outputPath,maskDiskDecr);
end
MakePath(outputPath);
% Parameters for a reconstruction algorithm using the GPU
cfg.ReconstructionDataId = rec_id;
cfg.ProjectionDataId = sino_id;
cfg.Options.GPUindex = 1; % Use GPU #1 for the reconstruction. (The default is #0.)
% Create the algorithm object from the configuration structure
alg_id = astra_mex_algorithm('create', cfg);

roiCircle = meshgrid((1:dimHorVol)-dimHorVol/2-0.5,(1:dimHorVol)-dimHorVol/2-0.5);
roiCircle = single(sqrt(roiCircle.^2+roiCircle'.^2) < 0.85*dimHorVol/2);
% Preallocation
rec = zeros([size(fbp.red),maxIter],'single');
error.fbp.full = zeros(1, maxIter);
error.fbp.red = zeros(1, maxIter);
error.res = zeros(1, maxIter);
% Iterate the algorithm one at a time, keeping track of errors
ecount = 0;
% Output folder name
outputPathIter = MakePath('%siterations_%s',outputPath,recType);
for nn = 1:maxIter
    PrintNum(nn)
    % Run a single iteration
    astra_mex_algorithm('iterate', alg_id, 1);
    % Get and save result
    rec(:,:,nn) = astra_mex_data2d('get_single', rec_id);
    %astra_mex_data2d('set',rec_id,FilterMedian(rec(:,:,nn)));
    filename = sprintf('%s%s_iter%04u',outputPathIter,recType,nn);
    if nn < 11 || mod(nn,2)
        WriteImage(filename,squeeze(rec(:,:,nn)),'tif')
    end
    % Error and residual
    error.fbp.full(nn) = sqrt(sumsqr((rec(:,:,nn) - fbp.full).*roiCircle));
    error.fbp.red(nn)  = sqrt(sumsqr((rec(:,:,nn) - fbp.red).*roiCircle));
    error.res(nn) = astra_mex_algorithm('get_res_norm', alg_id);
    % Abort loop if Nan or Inf occur
    if sum(sum(isnan(rec(:,:,nn))))>0 || sum(sum(isinf(rec(:,:,nn))))>0
        fprintf('\n Abort loop: NaN or Inf occured!')
                % Crop preallocated variables
                if nn < maxIter
                    rec(:,:,nn+1:end) = [];
                    error.fbp.full(nn+1:end) = [];
                    error.fbp.red(nn+1:end) = [];
                    error.res(nn+1:end) = [];
                end
                % Quit loop
                break
    end
    % Abort loop if error to full fbp increases again for more than 1/nn
    % maxIter
    if nn > minIter*(1-0.3)
        %if error.full.fbp.full(nn) > error.full.fbp.full(nn-2)
        if error.fbp.full(nn) > error.fbp.full(nn-2)
            ecount = ecount + 1;
            if ecount > 0.3*nn
                fprintf('\n Abort loop: Error ART-FBP_full increased since %u iterations!',ecount)
                % Crop preallocated variables
                if nn < maxIter
                    rec(:,:,nn+1:end) = [];
                    error.fbp.full(nn+1:end) = [];
                    error.fbp.red(nn+1:end) = [];
                    error.res(nn+1:end) = [];
                end
                % Quit loop
                break
            end
        else
            ecount = 0;
        end
    end
    
end
fprintf('\n Elapsed time: %g s\n',toc)

% Clean up. Note that GPU memory is tied up in the algorithm object,
% and main RAM in the data objects.
astra_mex_algorithm('delete', alg_id);
astra_mex_data2d('delete', rec_id);
astra_mex_data2d('delete', sino_id);
if maskDiskDecr > 0
    astra_mex_data2d('delete', mask_id);
end
%% Show plots, figures, images, ....%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if showFigs
    [roi.fbp, roi.art] = testAstraSaveImages(fbp,rec,error,outputPath,recType,showFigs);
    nimplay(rec(x,x,1:1:end),0)
else
    testAstraSaveImages(fbp,rec,error,outputPath,recType,showFigs);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n',hyphens)