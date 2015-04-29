function [recon_id, recon] = test_astra_fbp_cuda(proj_geom, vol_geom, sinogram, filter_type)

% Create a GPU based FBP reconstruction.
%
% proj_geom: projection geometry struct
% vol_geom: volume geometry struct
% sinogram: sinogram data OR sinogram identifier
% iterations: number of iterations to perform
% use_mask: use a reconstrucionmask? 'yes' or 'no'
% mask: mask data OR mask identifier.
% use_minc: use a minimum constraint? 'yes' or 'no'
% minc: minimum constraint value
% use_maxc: use a maximum constraint? 'yes' or 'no'
% maxc: maximum constraint value
% recon_id: identifier of the reconstruction data object as it is now stored in the astra-library
% recon: MATLAB data version of the reconstruction


if numel(sinogram) == 1
	sinogram_id = sinogram;
else
	sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, sinogram);
end

% create reconstruction object
recon_id = astra_mex_data2d('create', '-vol', vol_geom, 0);

% configure
cfg = astra_struct('FBP_CUDA');
%cfg.ProjectionGeometry = proj_geom;
%cfg.ReconstructionGeometry = vol_geom;
cfg.ProjectionDataId = sinogram_id;
cfg.ReconstructionDataId = recon_id;
% Use GPU index 1 for the reconstruction. The default is 0.
cfg.option.GPUindex = 1;
if nargin < 4
    cfg.FilterType = 'Ram-Lak';
else
    cfg.FilterType = filter_type;
end


% Create the algorithm object from the configuration structure
alg_id = astra_mex_algorithm('create', cfg);

% iterate
astra_mex_algorithm('iterate', alg_id, 0);

% return object
recon = astra_mex_data2d('get_single', recon_id);

% garbage collection
astra_mex_algorithm('delete', alg_id);
if numel(sinogram) ~= 1
	astra_mex_data2d('delete', sinogram_id);
end


