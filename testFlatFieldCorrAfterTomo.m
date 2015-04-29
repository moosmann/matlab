% BULLSHIT

dataFolder = '/mnt/tomoraid-LSDF/users/moosmann/CWI_DATA/ESRF_MI1079_ID19_July2011_inlineTomo/data/Xenopus_4cell_20keV/';

%% FLAT
ref = double(pmedfread([dataFolder 'refHST1599.edf'])');
%% tomography
% parameters
RotAxisPos = 1067;
NumProjTomo = 1599;
angles = 2*pi/NumProjTomo * ((1:NumProjTomo)-1);
% reconstruction
im = repmat(squeeze(ref(1024,:)'),[1 NumProjTomo]);
sinoFlat = (RotAxisSymmetricCropping(im(:,1:NumProjTomo)',RotAxisPos));
volFlat = astra_make_reco(sinoFlat,angles,'FBP_CUDA',1);

%% SAMPLE
im = MakeSino(1024,dataFolder,'Xenopus1_20keV_*',[],[],1:1599,0.01)';
%% tomography
% parameters
NumProjTomo = 1599;
angles = 2*pi/NumProjTomo * ((1:NumProjTomo)-1);
% reconstruction
sino = (RotAxisSymmetricCropping(im(:,1:NumProjTomo)',RotAxisPos));
vol = astra_make_reco(sino,angles,'FBP_CUDA',1);