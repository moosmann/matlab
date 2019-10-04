function vol = astra_parallel3D( par, sino)
% Parallel backprojection of 2D or 3D sinograms using ASTRA's
% parallel 3D geometry with vector notation. 
%
% sino: 2D-or-3D array.
% angles: scalar or vector. Default: pi. If calar it is the angular range
%   covered during one tomogram and the angles are computed as angles * (0:num_proj-1) /
%   num_proj. If vector it is the angles of the projections. If scalar the 
% rotation_axis_offset: scalar. Default: 0. Offset to the position of the
%   rotation axis. The rotation axis position is assumed to be the detector
%   center, size(sino,1)/2, shifted  by the rotation axis offset.
% vol_shape: shape of the reconstructed volume. Default: horizontal and
%   vertical number of voxel is given by the number of pixels of sinogram
%   along the first and second direction, respectively.
% vol_size: size of the reconstructed volume
% pixel_size: scalar or 2-component vector. Default: 1. Size of a detector
%   pixel. If scalar square pixels are assumed.
% link_data: boolean. Default: 0. If 0 ASTRA and MATLAB use their own
%   memory. If 1 ASTRA's data objects are references to MATLAB arrays.
%   Changes by ASTRA are visible to MATLAB. Changes by MATLAB creates a copy
%   of the data object and are not visible to the data object. Take care if 
%   using data links.
% tilt: scalar, tilt of rotation axis perpendicular to the beam. this accounts for
%   a rotation of the camera
% gpu_index: scalar, MATLAB index of GPU device to use. default: [], uses
%   all available GPUs. Matlab index notation starts from 1, ASTRA index
%   starts from 0. Here, MATLAB index notation is used.
% tilt_lamino : scalar, default: 0. tilt of rotation axis towards forward
%   beam direction in radians.
%
% For GPUs the only interpolation method available in ASTRA is the Josehp
% kernel.
%
% Written by Julian Moosmann
% First version: 2016-10-5. Last modification: 2017-10-30

%% TODO: test double precision support
%% TODO: proper Ram-Lak filter for tilted axis
%% TODO: check normalization factor pi / (2 * # angles)

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
angles = assign_from_struct( par, 'angles', pi );
rotation_axis_offset = assign_from_struct( par.rot_axis, 'offset', 0);
angle_offset = assign_from_struct( par.rot_angle, 'offset', 0 );
vol_shape = assign_from_struct( par, 'vol_shape', [size( sino, 1), size( sino, 1), size(sino, 3) ] );
vol_size = assign_from_struct( par, 'vol_size', [] );
pixel_size = assign_from_struct( par, 'astra_pixel_size', [1 1] );
tilt = assign_from_struct( par, 'tilt_camera', 0 );
tilt_lamino = assign_from_struct( par, 'tilt_lamino', 0 );
link_data = assign_from_struct( par, 'astra_link_data', 0 );
gpu_index = assign_from_struct( par, 'astra_gpu_index', [] );
algorithm = assign_from_struct( par, 'algorithm', 'fbp' );
iterations = assign_from_struct( par, 'iterations', 100);
MinConstraint = assign_from_struct( par.sirt, 'MinConstraint', [] );
MaxConstraint = assign_from_struct( par.sirt, 'MaxConstraint', [] );
vert_shift = assign_from_struct( par, 'vert_shift', [] );
 
%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

angles = double( angles );
%rotation_axis_offset = double( rotation_axis_offset );

% GPU
if isempty( gpu_index )
    astra_mex('set_gpu_index', 0:gpuDeviceCount - 1);
else
    astra_mex('set_gpu_index', gpu_index - 1);
end

%% Detector geometry
det_col_count = size( sino, 1);
det_row_count = size( sino, 3);
num_proj = size( sino, 2);
if numel(angles) == 1
    angles = angles * (0:num_proj-1) / num_proj;
end
if numel( pixel_size ) == 1
    pixel_size(2) = pixel_size;
end    
DetectorSpacingX = pixel_size(1);
DetectorSpacingY = pixel_size(2);
if numel( angles ) ~= size( sino, 2)
    error('Size of ANGLES and size of sinogram do not match.')
end

% Create geometry vector
vectors = zeros( numel(angles), 12);
for nn = 1:num_proj
    
    theta = angles( nn ) + angle_offset;
    
    % Lateral shift
    if isequal( numel( rotation_axis_offset ), 1 )
        rao = rotation_axis_offset;
    else
        rao = rotation_axis_offset(nn);
    end    

    % vertical shift
    if isempty( vert_shift )
        z = 0;
    else
        z = vert_shift(nn);
    end

    % source / ray direction
    %% CHECK
    vectors(nn,1) =  sin( theta );
    vectors(nn,2) = -cos( theta );
    vectors(nn,3) = z - sin(tilt_lamino);

    % center of detector
    vectors(nn,4) = -rao * cos( theta );
    vectors(nn,5) = -rao * sin( theta );
    vectors(nn,6) = z + sin( tilt_lamino );

    % vector from detector pixel (0,0) to (0,1)
    vectors(nn,7) = cos( tilt ) * cos( theta ) * DetectorSpacingX;
    vectors(nn,8) = cos( tilt) * sin( theta ) * DetectorSpacingX;
    vectors(nn,9) = cos( tilt_lamino ) * sin( tilt ) * DetectorSpacingX;

    % vector from detector pixel (0,0) to (1,0)
    vectors(nn,10) = -sin( tilt) * cos( theta ) * DetectorSpacingY;
    vectors(nn,11) = -sin( tilt) * sin( theta ) * DetectorSpacingY;
    vectors(nn,12) = cos( tilt_lamino) * cos(tilt) * DetectorSpacingY;

end

%% Projection geometry
proj_geom = astra_create_proj_geom('parallel3d_vec',  det_row_count, det_col_count, vectors);

%% Volume geometry
% Volume shape: y, x, z
%vol_size_astra = [vol_shape(2), vol_shape(1), vol_shape(3)];
row_count = vol_shape(2);
col_count = vol_shape(1);
slice_count = vol_shape(3);
% Volume size
if isempty( vol_size )
    vol_size = [-vol_shape(1)/2, vol_shape(1)/2, ...
        -vol_shape(2)/2, vol_shape(2)/2, -vol_shape(3)/2, vol_shape(3)/2];
end 
min_x = vol_size(3);
max_x = vol_size(4);
min_y = vol_size(1);
max_y = vol_size(2);
min_z = vol_size(5);
max_z = vol_size(6);
% Volume geometry object
vol_geom = astra_create_vol_geom(row_count, col_count, slice_count, min_x, max_x, min_y, max_y, min_z, max_z);

% Normalize sino instead of volume
sino = pi / 2 / length(angles) * sino;

%% Sinogram object
if link_data
    sino_id = astra_mex_data3d_c('link', '-proj3d', proj_geom, sino, 1, 0);
else
    sino_id = astra_mex_data3d('create', '-proj3d', proj_geom, sino);
end

%% Volume object
if link_data
    vol = zeros(vol_shape, 'single');
    vol_id = astra_mex_data3d_c('link', '-vol', vol_geom, vol, 1, 0);    
else
    vol_id = astra_mex_data3d('create', '-vol', vol_geom);
end

%% ASTRA config struct
switch lower( algorithm )
    case 'fbp'
        cfg = astra_struct('BP3D_CUDA');
    case 'sirt'
        cfg = astra_struct('SIRT3D_CUDA');
        cfg.option.MinConstraint = MinConstraint;
        cfg.option.MaxConstraint = MaxConstraint;
    case 'cgls'
        cfg = astra_struct('CGLS3D_CUDA');        
end
cfg.ProjectionDataId = sino_id;
cfg.ReconstructionDataId = vol_id;

%% ASTRA create algorithm object from configuration struct
bp_id = astra_mex_algorithm('create', cfg);

%% Backprojection: iterate algorithm 
if strcmpi( algorithm, 'fbp' )
    astra_mex_algorithm('iterate', bp_id, 1);
else
    astra_mex_algorithm('iterate', bp_id, iterations);
end
astra_mex_algorithm('delete', bp_id);
astra_mex_data3d('delete', sino_id)

%% Fetch data from ASTRA memory
if ~link_data
    vol = astra_mex_data3d('get_single', vol_id);
end
astra_mex_data3d('delete', vol_id)
astra_mex_projector3d('clear')

% Required for adjoint?
%vol = pi / 2 / numel(angles) * vol;

end
