outPath = '/ufs/moosmann/Pictures';

%% Halo
if 0
    dataPath = '/export/scratch1/moosmann/ESRF_MI1079_ID19_July2011_inlineTomo/vol/Xenopus_4cell_20keV/sino_intProjInc1__1600x1962x0100';
    sino = imread(sprintf('%s/slice_0050.tif',dataPath));
    
    recHalo = astra_make_reco(sino,2*pi*(0:1599)/1600);
    recNoHalo = astra_make_reco(SubtractMean(sino),2*pi*(0:1599)/1600);
    
    subOutPath = MakePath([outPath '/halo']);
    h = @(im,imString) WriteImage(sprintf('%sxeno4cell_rec_%s.png',subOutPath,imString),im);
    h(recHalo,'sinoNoMeanSub')
    h(recNoHalo,'sinoMeanSub')
end

%% CGLS, SIRT vs FBP
x = 500:1500;y = 500:1700;
dataPath2 = '/export/scratch1/moosmann/APS_32ID-C_InVivoImaging_GUP34879_2013-07-30/vol/Xenopus_inVivo/Jul29_15-10_urea_stage27p0_30p0keV_0700mm_15ms_0500proj_scantime20s_deadtime8min/tomo01/';
subOutPath = MakePath([outPath '/xeno27_aps_32id_fpb_nnfpb_cgls_sirt']);
hroi = @(im) im(x,y);

% read data
% ProjInc 1
intInc1_fbp  = imread([dataPath2 'tomo_intProjInc1_fbp__1966x1966x0100/slice_0050.tif']);
intInc1_cgls = imread([dataPath2 'tomo_intProjInc1_CGLS__1966x1966x0100/slice_0050.tif']);
intInc1_sirt = imread([dataPath2 'tomo_intProjInc1_SIRT__1966x1966x0100/slice_0050.tif']);
intInc1_fbp_tieRP25 = imread([dataPath2 'tomo_intProjInc1_fbp_tieRP25__1966x1966x0100/slice_0050.tif']);
intInc1_cgls_tieRP25 = imread([dataPath2 'tomo_intProjInc1_CGLS_tieRP25__1966x1966x0100/slice_0050.tif']);
intInc1_sirt_tieRP25 = imread([dataPath2 'tomo_intProjInc1_SIRT_tieRP25__1966x1966x0100/slice_0050.tif']);
tieRP25Inc1_fbp = imread([dataPath2 'tomo_tieRP25ProjInc1_fbp__1966x1966x0100/slice_0050.tif']);
% ProjInc 2
intInc2_fbp   = imread([dataPath2 'tomo_intProjInc2_fbp__1966x1966x0100/slice_0050.tif']);
intInc2_cgls  = imread([dataPath2 'tomo_intProjInc2_CGLS__1966x1966x0100/slice_0050.tif']);
intInc2_sirt  = imread([dataPath2 'tomo_intProjInc2_SIRT__1966x1966x0100/slice_0050.tif']);
intInc2_nnfbp = imread([dataPath2 'tomo_intProjInc2_nnfbp__1966x1966x0100/slice_0050.tif']);
intInc2_fbp_tieRP25   = imread([dataPath2 'tomo_intProjInc2_fbp_tieRP25__1966x1966x0100/slice_0050.tif']);
intInc2_cgls_tieRP25  = imread([dataPath2 'tomo_intProjInc2_CGLS_tieRP25__1966x1966x0100/slice_0050.tif']);
intInc2_sirt_tieRP25  = imread([dataPath2 'tomo_intProjInc2_SIRT_tieRP25__1966x1966x0100/slice_0050.tif']);
intInc2_nnfbp_tieRP25 = imread([dataPath2 'tomo_intProjInc2_nnfbp_tieRP25__1966x1966x0100/slice_0050.tif']);
% ProjInc 3
intInc3_fbp         = imread([dataPath2 'tomo_intProjInc3_fbp__1966x1966x0100/slice_0050.tif']);
intInc3_fbp_tieRP25 = imread([dataPath2 'tomo_intProjInc3_fbp_tieRP25__1966x1966x0100/slice_0050.tif']);

% Write images
% ProjInc1
% Int, full
h = @(im,imInd,imString) WriteImage(sprintf('%s/%02u_intProjInc1_full_%s',subOutPath,0+imInd,imString),im(x,y),'png');
h(intInc1_fbp,   1,'fbp')
h(intInc1_sirt,  2,'sirt')
h(intInc1_cgls,  3,'cgls')
% Phase, full
h = @(im,imInd,imString) WriteImage(sprintf('%s/%02u_intProjInc1_tieRP25_full_%s',subOutPath,10+imInd,imString),im(x,y),'png');
h(intInc1_fbp_tieRP25,   1,'fbp')
h(intInc1_sirt_tieRP25,  2,'sirt')
h(intInc1_cgls_tieRP25,  3,'cgls')
% ProjInc2
% Int, full
h = @(im,imInd,imString) WriteImage(sprintf('%s/%02u_intProjInc2_full_%s',subOutPath,20+imInd,imString),im(x,y),'png');
h(intInc2_fbp,   1,'fbp')
h(intInc2_sirt,  2,'sirt')
h(intInc2_cgls,  3,'cgls')
h(intInc2_nnfbp, 4,'nnfbp')
% Phase, full
h = @(im,imInd,imString) WriteImage(sprintf('%s/%02u_intProjInc2_tieRP25_full_%s',subOutPath,30+imInd,imString),im(x,y),'png');
h(intInc2_fbp_tieRP25,  1,'fbp')
h(intInc2_sirt_tieRP25, 2,'sirt')
h(intInc2_cgls_tieRP25, 3,'cgls')
h(intInc2_nnfbp_tieRP25,4,'nnfbp')
% Phase, roi
x2 = 730+(1:256); y2 = 730+(1:256);
h = @(im,imInd,imString) WriteImage(sprintf('%s/%02u_intProjInc2_tieRP25_roi_%s',subOutPath,40+imInd,imString),im(x2,y2),'png');
h(intInc2_fbp_tieRP25,  1,'fbp')
h(intInc2_sirt_tieRP25, 2,'sirt')
h(intInc2_cgls_tieRP25, 3,'cgls')
h(intInc2_nnfbp_tieRP25,4,'nnfbp')

%% Phase tomo VS tomo phase
if 0
    % read data
    % xeno 4cell
    dataPath0 = '/export/scratch1/moosmann/ESRF_MI1079_ID19_July2011_inlineTomo/vol/Xenopus_4cell_20keV/';
    qp25bf01 = imread([dataPath0 'tomo_qpRP25BF01Proj1__1962x1962x0100/slice_0050.tif']);
    fbp_qp25bf01 = imread([dataPath0 'tomo_intProjInc1_fbp_qpRP25BF01__1962x1962x0100/slice_0050.tif']);
    %xeno 27
    dataPath2 = '/export/scratch1/moosmann/APS_32ID-C_InVivoImaging_GUP34879_2013-07-30/vol/Xenopus_inVivo/Jul29_15-10_urea_stage27p0_30p0keV_0700mm_15ms_0500proj_scantime20s_deadtime8min/tomo01/';
    tieRP25Inc1_fbp = imread([dataPath2 'tomo_tieRP25ProjInc1_fbp__1966x1966x0100/slice_0050.tif']);
    intInc1_fbp_tieRP25 = imread([dataPath2 'tomo_intProjInc1_fbp_tieRP25__1966x1966x0100/slice_0050.tif']);
    
    % Write images
    x = 500:1500;y = 500:1700;
    subOutPath = MakePath([outPath '/phaseTomo_vs_tomoPhase']);
    % full
    h = @(im,imInd,imString) WriteImage(sprintf('%s/%02u_xeno27_projInc1_full_%s',subOutPath,0+imInd,imString),im,'png');
    h(tieRP25Inc1_fbp,   1,'tieRP25_fbp')
    h(intInc1_fbp_tieRP25,  2,'fbp_tieRP25')
    % roi
    h = @(im,imInd,imString) WriteImage(sprintf('%s/%02u_xeno27_projInc1_roi_%s',subOutPath,10+imInd,imString),im(x,y),'png');
    h(tieRP25Inc1_fbp,   1,'tieRP25_fbp')
    h(intInc1_fbp_tieRP25,  2,'fbp_tieRP25')
    
    % full
    h = @(im,imInd,imString) WriteImage(sprintf('%s/%02u_xeno4cell_projInc1_full_%s',subOutPath,20+imInd,imString),im,'png');
    h(qp25bf01,      1,'qpRP25BF01_fbp')
    h(fbp_qp25bf01,  2,'fbp_qpRP25BF01')
    % roi
    h = @(im,imInd,imString) WriteImage(sprintf('%s/%02u_xeno4cell_projInc1_roi_%s',subOutPath,30+imInd,imString),im(x,x),'png');
    h(qp25bf01,      1,'qpRP25BF01_fbp')
    h(fbp_qp25bf01,  2,'fbp_qpRP25BF01')
end