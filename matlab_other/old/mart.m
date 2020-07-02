
clear all
aclear

% Read sinogramm
sino = SubtractMean(getSino('xeno4cell',12));
% Parameter
[NumProjFull, dimHor] = size(sino);
dimHorVol = 1*dimHor;
angularIncrementInRad = 2*pi/1599;
NumProj = NumProjFull;
%angles = @(NumProj,angleInc) angleInc*angularIncrementInRad*((1:NumProj)-1);
angles =  angularIncrementInRad*((1:NumProj)-1);
% Creat volume and projection geometry
vol_geom = astra_create_vol_geom([dimHorVol,dimHorVol]);
projUsed = 1:NumProj;
proj_geom = astra_create_proj_geom('parallel',1.0,dimHor,angles(projUsed));
% Create the sinogram data object
%sino_id = astra_mex_data2d('create', '-sino', proj_geom, sino(projUsed,:));
% Create reconstruction
%[fbp_id, fbp] = test_astra_fbp_cuda(proj_geom,vol_geom,sino_id);
%astra_mex_data2d('delete', fbp_id);
bp_sino = astra_create_backprojection_cuda(sino,proj_geom,vol_geom);
s = sino(projUsed,:);
k = abs(FrequencyVector(dimHor));
s = ifft(repmat(k,[size(s,1) 1]).*fft(s,[],2),[],2);
% astra backprojection of filtered sinogram 's'
tic
bp = astra_create_backprojection_cuda(s,proj_geom,vol_geom);
toc

% iradon filtered backprojection of sinogram 'sino'
sino_ir = iradon(sino',360/1599*(0:1599),'linear','ram-lak',1,dimHor)*dimHor/2;
% iradon unfiltered backprojection of sinogram 'sino'
sino_ir_nof = iradon(sino',360/1599*(0:1599),'linear','none',1,dimHor)*dimHor/2;
% iradon unfiltered backprojection of filtered sinogram 's'
tic
s_ir = iradon(s',360/1599*(0:1599),'linear','none',1,dimHor)*dimHor/2;

toc
% projection of astra backprojection
[fp_id, fp] = astra_create_sino_cuda(SubtractMean(bp), proj_geom, vol_geom);

%astra_mex_data2d('delete', sino_id);
