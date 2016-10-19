function vol = astra_parallel3D(sino, angles, rotation_axis_offset, vol_shape, pixel_size, link_data)
% Parallel backprojection of 2D or 3D sinograms using ASTRA's
% parallel 3D geometry with vector notation. 
%
% sino: 2D-or-3D array.
% angles: scalar or vector. Default: pi. Angular range or array of angles
% of the projection. If scalar the angles are angles * (0:num_proj-1) /
% num_proj. 
% rotation_axis_offset: scalar. Default: 0. Offset to the position of the
% rotation axis. The rotation axis position is assumed to be the detector
% center, size(sino,1)/2, shifted  by the rotation axis offset.
% vol_shape: size of the to be reconstructed volume. Default: horizontal and
% vertical number of voxel is given by the number of pixels of sinogram
% along the first and second direction, respectively.
% pixel_size: 2-component vector. Default: [1 1]. Length a detector pixel.
% link_data: boolean. Default: 0. If 0 ASTRA and MATLAB use their own
% memory. If 1 ASTRA's data objects are references to MATLAB arrays.
% Changes inside by ASTRA are visible to MATLAB. Changes by MATLAB creates
% a copy of the data object and are not visible to the data object. Take if
% using data links.
%
% Written by Julian Moosmann
% First version: 2016-10-5. Last modification: 2016-10-12

%% TODO: test double precision support

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    angles = pi;
end
if nargin < 3
    rotation_axis_offset = 0;
end
if nargin < 4
    vol_shape = [size( sino, 1), size( sino, 1), size(sino, 3) ];
end
if nargin < 5
    pixel_size = [1, 1];
end
if nargin < 6
    link_data = 0;
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Detector geometry
det_col_count = size( sino, 1);
det_row_count = size( sino, 3);
num_proj = size( sino, 2);
if numel(angles) == 1
    angles = angles * (0:num_proj-1) / num_proj;
end
DetectorSpacingX = pixel_size(1);
DetectorSpacingY = pixel_size(2);
if numel( angles ) ~= size( sino, 2)
    error('Size of ANGLES and size of sinogram do not match.')
end

% Create geometry vector
vectors = zeros( numel(angles), 12);
for nn = 1:num_proj
    
    theta = angles( nn );

    % source / ray direction
    vectors(nn,1) = sin( theta );
    vectors(nn,2) = -cos( theta );
    vectors(nn,3) = 0;

    % center of detector
    vectors(nn,4) = -rotation_axis_offset * cos( theta );
    vectors(nn,5) = -rotation_axis_offset * sin( theta );
    vectors(nn,6) = 0;

    % vector from detector pixel (0,0) to (0,1)
    vectors(nn,7) = cos( theta ) * DetectorSpacingX;
    vectors(nn,8) = sin( theta ) * DetectorSpacingX;
    vectors(nn,9) = 0;

    % vector from detector pixel (0,0) to (1,0)
    vectors(nn,10) = 0;
    vectors(nn,11) = 0;
    vectors(nn,12) = DetectorSpacingY;

end

% Projection geometry
proj_geom = astra_create_proj_geom('parallel3d_vec',  det_row_count, det_col_count, vectors);

% Volume geometry: y, x, z
vol_size_astra = [vol_shape(2), vol_shape(1), vol_shape(3)];
vol_geom = astra_create_vol_geom(vol_size_astra);

% Normalize sino instead of volume
sino = pi / 2 / length(angles) * sino;

% Sinogram object
if link_data
    sino_id = astra_mex_data3d_c('link', '-proj3d', proj_geom, sino, 1, 0);
else
    sino_id = astra_mex_data3d('create', '-proj3d', proj_geom, sino);
end

% Volume object
if link_data
    vol = zeros(vol_shape, 'single');
    vol_id = astra_mex_data3d_c('link', '-vol', vol_geom, vol, 1, 0);    
else
    vol_id = astra_mex_data3d('create', '-vol', vol_geom);
end

% ASTRA config struct
cfg = astra_struct('BP3D_CUDA');
cfg.ProjectionDataId = sino_id;
cfg.ReconstructionDataId = vol_id;

% ASTRA create algorithm object from configuration struct
bp_id = astra_mex_algorithm('create', cfg);

% Backprojection (one algorithm iteration)
astra_mex_algorithm('iterate', bp_id, 1);
astra_mex_algorithm('delete', bp_id);
astra_mex_data3d('delete', sino_id)

% Fetch data from ASTRA memory
if ~link_data
    vol = astra_mex_data3d('get_single', vol_id);
end
astra_mex_data3d('delete', vol_id)

% Normalize volume
%vol = pi / 2 / numel(angles) * vol;
