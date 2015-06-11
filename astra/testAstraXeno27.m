aclear
ParentPath = '/export/scratch1/moosmann/APS_32ID-C_InVivoImaging_GUP34879_2013-07-30/vol/Xenopus_inVivo/Jul29_15-10_urea_stage27p0_30p0keV_0700mm_15ms_0500proj_scantime20s_deadtime8min/tomo01';
%% Read sinogramm
%sino = SubtractMean(imread('slice_0050.tif'));
sino = Readstack(sprintf('%s/sino_intProjInc1__0499x1966x0100',ParentPath));

%% Parameter
projInc = 2;
[NumProjFull, dimHor] = size(sino(:,:,1));
dimHorVol = 1*dimHor;
angularIncrementInRad = pi/500;
NumProj = NumProjFull;
angles =  angularIncrementInRad*((1:NumProj)-1);

%% Creat volume and projection geometry, and sinogram
vol_geom = astra_create_vol_geom([dimHorVol,dimHorVol]);
projUsed = 1:projInc:NumProj;
proj_geom = astra_create_proj_geom('parallel',1.0,dimHor,angles(projUsed));

slices = size(sino,3):-1:1;
rectypes = {'CGLS','SIRT'};
for mm = 1:numel(rectypes)
    rectype = rectypes{mm};
    %% Create reconstruction
    for nn = numel(slices):-1:1
        PrintNum(nn);
        sliceNum = slices(nn);
        sino_id = astra_mex_data2d('create', '-sino', proj_geom, sino(projUsed,:,sliceNum));
        [rec_id, rec(:,:,nn)] = astra_create_reconstruction_cuda(sprintf('%s_CUDA',rectype),proj_geom,vol_geom,sino_id,200,0,0,0,0,0,0);
        astra_mex_data2d('delete', sino_id);
        astra_mex_data2d('delete', rec_id);
    end
    
    %% Save volume
    WriteVol(rec,sprintf('%s/tomo_intProjInc%u_%s_',ParentPath,projInc,rectype))
    
end