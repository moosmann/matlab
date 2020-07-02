aclear;
switch 2
    case 1
        %% Parameters
        DataPath = '/export/scratch1/moosmann/APS_32-ID-C_2012-10-10_Life_Cell_Imaging/phase/Oct11_00-30_wildtype_stage11p0_34p5keV_0700mm_15ms_0834proj_scantime50s_deadtime08min_20ms_open_40ms_close/3DstackProc_FiltSino_FDcor_tie_regPar2p50_noMeanSub/tomo00_tif';
        SinoPath = '/export/scratch1/moosmann/APS_32-ID-C_2012-10-10_Life_Cell_Imaging/phase/Oct11_00-30_wildtype_stage11p0_34p5keV_0700mm_15ms_0834proj_scantime50s_deadtime08min_20ms_open_40ms_close/3DstackProc_FiltSino_FDcor_tie_regPar2p50_noMeanSub/sino00_tif';
        % From folder name
        tExp = 15;
        tOpen = 20;
        tClose = 40;
        tCycle = tOpen + tClose;
        numProj = 834;
        tTomoFold = 50*1000;
        tTomo = (tCycle)*numProj;
        pixelsize = 1.1;
        distFromRA = 1000;
        %%
        angBlur = pi*tExp/tTomo;
        angIncFromParDeg = 0.21635;
        fprintf('Exposure time per image: %g ms\n',tExp)
        fprintf('Exposure time per tomograms: %g ms * %u projections = %g s\n',tExp,numProj,numProj*tExp/1000)
        fprintf('Open and close time: %g ms and %g ms\n',tOpen,tClose)
        fprintf('Time per tomograms = (open time + close time) * projections: %g s, From folder name: %g s\n',tTomo/1000,tTomoFold/1000)
        fprintf('Angular blurring: %g mrad = %g mdeg \n',angBlur*1000, angBlur*180/pi*1000)
        fprintf('Pixel size: %g micron\n',pixelsize)
        fprintf('Arc length blurring at %g micron afar from rot axis: %g micron\n',distFromRA,angBlur*distFromRA)
        fprintf('Angular increment: 180*/projections = %g degree, from par file: %g degree\n',180/(numProj-1),angIncFromParDeg)
        % Read sino
        sino = imread(sprintf('%s/slice_0700.tif',SinoPath));
    case 2
         %% Parameters
        % ANGLE_BETWEEN_PROJECTIONS = 0.36 # Increment angle in degrees
        % ROTATION_AXIS_POSITION = 988.500000 # Position in pixels
        SinoPath = '/export/scratch1/moosmann/APS_32ID-C_InVivoImaging_GUP34879_2013-07-30/int/Xenopus_inVivo/Jul28_14-50_urea_stage17p0_30p0keV_0400mm_20ms_0500proj_scantime20s_deadtime8min/sino02_tif';
        % From folder name
        tExp = 20;
        numProj = 499;
        tTomoFold = 20*1000;        
        pixelsize = 1.1;
        distFromRA = 800;
        %%
        angBlur = pi*tExp/tTomoFold;
        angIncFromParDeg = 0.36;
        fprintf('Exposure time per image: %g ms\n',tExp)
        fprintf('Exposure time per tomograms: %g ms * %u projections = %g s\n',tExp,numProj,numProj*tExp/1000)
        fprintf('Time per tomograms from folder name: %g s\n',tTomoFold/1000)
        fprintf('Angular blurring: %g mrad = %g mdeg \n',angBlur*1000, angBlur*180/pi*1000)
        fprintf('Pixel size: %g micron\n',pixelsize)
        fprintf('Arc length blurring at %g micron afar from rot axis: %g micron\n',distFromRA,angBlur*distFromRA)
        fprintf('Angular increment: 180*/%u = %g degree, from par file: %g degree\n',numProj,180/(numProj-1),angIncFromParDeg)
        % Read sino
        sino = imread(sprintf('%s/slice_1001.tif',SinoPath));
end
OutputPath = '~/Pictures/Rotation_flyScan_vs_stepwise/experiment';
CheckAndMakePath(OutputPath);
%x = 300:1200;
dx = 512;
x = 1200 + (1:dx);
y = 600 + (1:dx);
sino = SubtractMean(sino);
%% Continous FBP
[numProjSino, dimHorVol] = size(sino);
angles = pi/numProj*(0:numProjSino - 1);
recoType = 'FBP_CUDA';
iterations = 1;
M = 100;
rec = zeros(dimHorVol,dimHorVol,M+1);
for nn = 0:M
    % Reco at different angles
    rec(:,:,nn+1) = astra_make_reco(sino,angles+nn/M*angBlur,recoType,iterations);
end
% Average over recos
recs = sum(rec,3)/(M+1);
% Pick single 'discrete' reco
rec0 = rec(:,:,ceil((M+1)/2));
%% Save images
% full
dynRange = [ max([min2(rec0) min2(recs)]) min([max2(rec0) max2(recs)])];
hf = @(im,imInd,imString,dynRange) WriteImage(sprintf('%s/%02u_full_%s.png',OutputPath,imInd,imString),SetDynRange(im,dynRange));
hf(rec0,11,'recoSingle',dynRange)
hf(recs,12,sprintf('recoSuperposOf%03u',M+1),dynRange);
hf(abs(rec0-recs),13,'diffMap',[0 0.8*max(abs(rec0-recs)) ] );
% roi
dynRangeRoi = [ max([min2(rec0(x,y)) min2(recs(x,y))]) min([max2(rec0(x,y)) max2(recs(x,y))])];
hr = @(im,imInd,imString,dynRange) WriteImage(sprintf('%s/%02u_roi_%s.png',OutputPath,imInd,imString),SetDynRange(im(x,y),dynRange));
hr(rec0,21,'recoSingle',dynRangeRoi);
hr(recs,22,sprintf('recoSuperposOf%03u',M+1),dynRangeRoi);
hr(abs(rec0-recs),23,'diffMap',[0 0.8*max(abs(rec0(x,y)-recs(x,y))) ] );

%% Continuous iteration
M = 3;
sino2 = zeros(M*numProjSino,dimHorVol);
angles2 = zeros(1,M*numProjSino);

for nn = 1:M
    sino2(nn:M:end,:) = sino;
    if M > 1
        angles2(nn:M:end) = angles(:) + (nn-1)/(M-1)*angBlur;
    else
        angles2 = angles;
    end
end
[numProjSino2, dimHorVol] = size(sino2);

vol_geom = astra_create_vol_geom([dimHorVol,dimHorVol]);
proj_geom = astra_create_proj_geom('parallel',1.0,dimHorVol,angles2);

% Create the sinogram data object
sino_id = astra_mex_data2d('create', '-sino', proj_geom, sino2);

% Create a data object for the reconstruction
rec_id = astra_mex_data2d('create', '-vol', vol_geom,0);
cfg = astra_struct('CGLS_CUDA');
% Parameters for a reconstruction algorithm using the GPU
cfg.ReconstructionDataId = rec_id;
cfg.ProjectionDataId = sino_id;
cfg.Options.GPUindex = 2; % Use GPU #1 for the reconstruction. (The default is #0.)
% Create the algorithm object from the configuration structure
alg_id = astra_mex_algorithm('create', cfg);
numIter = 200;
rec = zeros(dimHorVol,dimHorVol,numIter);
for nn = 1:numIter
    % Run reco
    astra_mex_algorithm('iterate', alg_id, 1);
    % Get result
    rec(:,:,nn) = astra_mex_data2d('get_single', rec_id);
end


% Clean up
astra_mex_algorithm('delete', alg_id);
astra_mex_data2d('delete', rec_id);
astra_mex_data2d('delete', sino_id);

itool(fbp,x,6),itool(rec,x,numIter)