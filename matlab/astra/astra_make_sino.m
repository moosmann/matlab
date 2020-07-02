function [sino, angles, vol_geom, proj_geom, cfg] = astra_make_sino(data, numPixVol, numPixSino)
% Written for Matlab Script 'ContRotRecoSim'.
% 
% Written by Julian Moosmann, last version: 2013-12-10
%
% [sino, angles, vol_geom, proj_geom, cfg] = astra_make_sino(data, numPixVol, numPixSino)

%% Defaults %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    % Number of volume pixels
    numPixVol = max(size(data));
end
if nargin < 3
    % Number of projections
    M = 10^(abs(floor(log10(numPixVol))-1));
    numPixSino = M*ceil(numPixVol*pi/2/M);
end

%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Geometries
% Create volume geometry
vol_geom = astra_create_vol_geom(numPixVol, numPixVol);

% Create a parallel beam geometry with 180 angles between 0 and pi, and
% 384 detector pixels of width 1.
angles = linspace2(0,pi,numPixSino);
proj_geom = astra_create_proj_geom('parallel', 1.0, numPixVol, angles);

%% Create sinogram
% store volume
vol_id = astra_mex_data2d('create','-vol', vol_geom, data);

% store sino
sino_id = astra_mex_data2d('create','-sino', proj_geom, 0);
    
% create struct
cfg = astra_struct('FP_CUDA');
cfg.ProjectionDataId = sino_id;
cfg.VolumeDataId = vol_id;
cfg.option.GPUindex = 1;

% create sinogram
alg_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('iterate', alg_id);
sino = astra_mex_data2d('get',sino_id);

% Free memory
astra_mex_data2d('delete', vol_id);
astra_mex_algorithm('delete', alg_id);
astra_mex_data2d('delete', sino_id);