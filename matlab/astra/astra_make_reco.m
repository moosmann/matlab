function [vol, vol_geom, proj_geom] = astra_make_reco(sino,angles,recoType,iterations)
% 
% Written for Matlab Script 'ContRotRecoSim'.
% 
% Written by Julian Moosmann, last version: 2013-12-10
%
% [vol, vol_geom, proj_geom] = astra_make_reco(sino,angles,recoType,iterations)

if nargin < 3 
    recoType = 'FBP_CUDA';
end
if nargin < 4
    iterations = 150;
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmpi(recoType(1:3),'FBP')
    iterations = 1;
end

dimHorVol = size(sino,1);
% Create geometry
vol_geom = astra_create_vol_geom(dimHorVol, dimHorVol);
% Create geometry
proj_geom = astra_create_proj_geom('parallel',1.0,dimHorVol,angles);
% Create the sinogram data object
sino_id = astra_mex_data2d('create', '-sino', proj_geom, sino);
% Create a data object for the reconstruction
vol_id = astra_mex_data2d('create', '-vol', vol_geom,0);
% Create and configure struct
cfg = astra_struct(recoType);
% Parameters for a reconstruction algorithm using the GPU
    cfg.ProjectionDataId = sino_id;
cfg.ReconstructionDataId = vol_id;
cfg.option.GPUindex = 1; 
if strcmpi(recoType,'FBP')
    cfg.FilterType = 'Ram-Lak';
end
% Create the algorithm object from the configuration structure
alg_id = astra_mex_algorithm('create', cfg);
% Run reco
astra_mex_algorithm('iterate', alg_id, iterations);
% Get result
vol = astra_mex_data2d('get', vol_id);
% Clean up
astra_mex_algorithm('delete', alg_id);
astra_mex_data2d('delete', vol_id);
astra_mex_data2d('delete', sino_id);