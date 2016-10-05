function vol = astra_parallel3D(sino, angles, rotationAxisOffset, volSize, pixelSize)
% Parallel backprojection of 2D or 3D sinograms using ASTRA's
% parallel 3D geometry with vector notation. 
%
% sino: 2D-or-3D array.
% angles: scalar or vector. Default: pi. Angular range or array of angles
% of the projection. If scalar the angles are angles * (0:num_proj-1) /
% num_proj. 
% rotationAxisOffset: scalar. Default: 0. Offset to the position of the
% rotation axis. The rotation axis position is assumed to be the detector
% center, size(sino,1)/2, shifted  by the rotation axis offset.
% volSize: size of the to be reconstructed volume. Default: horizontal and
% vertical number of voxel is given by the number of pixels of sinogram
% along the first and second direction, respectively.
% pixelSize: 2-component vector. Default: [1 1]. Length a detector pixel.
%
% Written by Julian Moosmann
% First version: 2016-10-5. Last modification: 2016-10-05

%% TODO: add support rotation axis offset

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    angles = pi;
end
if nargin < 3
    rotationAxisOffset = 0;
end
if nargin < 4
    volSize = [size( sino, 1), size( sino, 1), size(sino, 3) ];
end
if nargin < 5
    pixelSize = [1, 1];
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Detector geometry
det_col_count = size( sino, 1);
det_row_count = size( sino, 3);
num_proj = size( sino, 2);
if numel(angles) == 1
    angles = angles * (0:num_proj-1) / num_proj;
end
DetectorSpacingX = pixelSize(1);
DetectorSpacingY = pixelSize(2);
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
    vectors(nn,4) = -rotationAxisOffset * cos( theta );
    vectors(nn,5) = -rotationAxisOffset * sin( theta );
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
vol_geom = astra_create_vol_geom( volSize(2), volSize(1), volSize(3) );

% Volume object
vol_id = astra_mex_data3d('create', '-vol', vol_geom);

% Sinogram object
sino_id = astra_mex_data3d('create', '-proj3d', proj_geom, sino);

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
vol = astra_mex_data3d('get_single', vol_id);
astra_mex_data3d('delete', vol_id)

% Normalize volume
vol = pi / 2 / length(angles) * vol;
