
%%
parpath = '/asap3/petra3/gpfs/p05/2017/data/11003440/processed/syn20_62R_Mg10Gd_12w/reco/pd/toforward/';
preco = [ parpath 'cropped_binned' ];
pseg = [ parpath 'segmented' ];
voxel = [4.8051 4.8051 4.8051];

fprintf( '\n voxel size : [%f %f %f]', voxel )

%% 
parpath = '/asap3/petra3/gpfs/p05/2017/data/11003440/processed/syn32_99R_Mg10Gd_4w/reco/toforward/';
preco = 'syn32_original';
pseg = 'syn32_semented';
voxel = [0.0031496 0.0031496 0.00629922];

parpath = '/asap3/petra3/gpfs/p05/2017/data/11003440/processed/syn33_80R_Mg10Gd_8w/reco/toforward';
preco = 'syn33_original';
pseg = 'syn33_semented';
voxel = [0.0015748 0.0015748 0.00314961];

parpath = '/asap3/petra3/gpfs/p05/2017/data/11003440/processed/syn39_75L_Mg5Gd_8w/reco/toforward';
preco = 'syn39_original';
pseg = 'syn39_semented';
voxel = [0.0015748 0.0015748 0.00314961];
