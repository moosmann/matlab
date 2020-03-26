function vol = astra_parallel2D( tomo, sino, gpu_index, rotation_axis_offset)
% Slicewise parallel backprojection of 2D or 3D sinograms using ASTRA.
%
% ARGUMENTS
% tomo : parameter struct with fields:
%   angles: scalar or vector. Default: pi. If calar it is the angular range
%       covered during one tomogram and the angles are computed as angles * (0:num_proj-1) /
%       num_proj. If vector it is the angles of the projections. If scalar the 
%   rotation_axis_offset: scalar. Default: 0. Offset to the position of the
%       rotation axis. The rotation axis position is assumed to be the detector
%       center, size(sino,1)/2, shifted  by the rotation axis offset.
%   vert_shift : vector of vertical shifts for sprial CT, default: []
%   angle_offet : scalar, global angular offset to the reconstruction
%   vol_shape: shape of  reconstructed volume. Default: horizontal and
%       vertical number of voxel is given by the number of pixels of sinogram
%       along the first and second direction, respectively.
%   vol_size: size of reconstructed volume, default: [], then inferred from
%       vol_shape.
%   pixel_size: scalar or 2-component vector. Default: 1. Size of a detector
%       pixel. If scalar square pixels are assumed.
%   tilt: scalar, tilt of rotation axis perpendicular to the beam. this accounts for
%       a rotation of the camera
%   tilt_lamino : scalar, default: 0. tilt of rotation axis towards forward
%       beam direction in radians.
%   algorithm : string, see p05_reco
%   iterations : scalar, number of iteration for iterative methods,
%      default: 100
%   MinConstraint : scalar, for 
%   MaxConstraint : scalar, for    
%   gpu_index: scalar, MATLAB index of GPU device to use. default: [], uses
%       all available GPUs. Matlab index notation starts from 1, ASTRA index
%       starts from 0. Here, MATLAB index notation is used.
%   link_data: boolean. Default: 0. If 0 ASTRA and MATLAB use their own
%       memory. If 1 ASTRA's data objects are references to MATLAB arrays.
%       Changes by ASTRA are visible to MATLAB. Changes by MATLAB creates a copy
%       of the data object and are not visible to the data object. Take care if 
%       using data links.
% sino: 2D-or-3D array.
%
% For GPUs the only interpolation method available in ASTRA is the Josehp
% kernel.
%
%% Note: For usage with GPU code, the volume must be centered around the
%% origin and pixels must be square. This is not always explicitly checked
%% in all functions, so not following these requirements may have
%% unpredictable results.
%
% Written by Julian Moosmann

if nargin < 3
    gpu_index = [];
end
if nargin < 4
    rotation_axis_offset = [];
end

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
angles = assign_from_struct( tomo, 'angles', pi );
vol_shape = assign_from_struct( tomo, 'vol_shape', [size( sino, 1), size( sino, 1), size(sino, 3) ] );
vol_size = assign_from_struct( tomo, 'vol_size', [] );
pixel_size = assign_from_struct( tomo, 'astra_pixel_size', 1 );
tilt = assign_from_struct( tomo, 'rot_axis_tilt_camera', 0 );
algorithm = assign_from_struct( tomo, 'algorithm', 'fbp' );
iterations = assign_from_struct( tomo, 'iterations', 100);
scan_position = assign_from_struct( tomo, 'scan_position', 0);
angle_offset = assign_from_struct( tomo, 'rot_angle_offset', 0 );
MinConstraint = assign_from_struct( tomo, 'MinConstraint', [] );
MaxConstraint = assign_from_struct( tomo, 'MaxConstraint', [] );
%vert_shift = assign_from_struct( tomo, 'vert_shift', [] );
astra_gpu_index = assign_from_struct( tomo, 'astra_gpu_index', [] );
if isempty( rotation_axis_offset )
    rotation_axis_offset = assign_from_struct( tomo, 'rot_axis_offset', 0);
end

%% TODO: Spiral CT using interpolation

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if tilt ~= 0
    error( 'Rotation axis tilt is not supported for 2D slicewise reconstructions' )
end

angles = double( angles );
%rotation_axis_offset = double( rotation_axis_offset );

% GPU
if isempty( gpu_index )
    if ~isempty( astra_gpu_index )
        astra_mex('set_gpu_index', astra_gpu_index - 1);
    else
        astra_mex('set_gpu_index', 0:gpuDeviceCount - 1);
    end
else
    nn = astra_gpu_index( 1 + mod( gpu_index, numel( astra_gpu_index ) ) );
    astra_mex('set_gpu_index', nn - 1);
end

%% Detector geometry
det_col_count = size( sino, 2);
num_proj = size( sino, 1);
if numel(angles) == 1
    angles = angles * (0:num_proj-1) / num_proj;
end
if numel( pixel_size ) == 2
    error( 'astra_pixel_size (%f,%f) does not fit geometry', pixel_size )
end    
DetectorSpacingX = pixel_size(1);
if numel( angles ) ~= size( sino, 1)
    error('Size of ANGLES and size of sinogram do not match.')
end

% Create geometry vector
vectors = zeros( numel(angles), 6);
for nn = 1:num_proj
    
    theta = angles( nn ) + angle_offset;
    
    % lateral shift
    if isequal( numel( rotation_axis_offset ), 1 )
        rao = rotation_axis_offset;
    else
        rao = rotation_axis_offset(nn);
    end

    % Scan position
    if ~isscalar( scan_position )
        rao = rao + scan_position(nn);
    end
    
    % source / ray direction
    %% CHECK
    vectors(nn,1) = + sin( theta );
    vectors(nn,2) = - cos( theta );

    % center of detector
    vectors(nn,3) = -rao * cos( theta );
    vectors(nn,4) = -rao * sin( theta );

    % vector from detector pixel 0 to 1
    vectors(nn,5) = cos( theta ) * DetectorSpacingX;
    vectors(nn,6) = sin( theta ) * DetectorSpacingX;
end

%% Projection geometry
proj_geom = astra_create_proj_geom('parallel_vec', det_col_count, vectors);

%% Volume geometry
% Volume shape: y, x
row_count = vol_shape(2);
col_count = vol_shape(1);
% Volume size
if isempty( vol_size )
    vol_size = [-vol_shape(1)/2, vol_shape(1)/2, -vol_shape(2)/2, vol_shape(2)/2];
end 
min_x = vol_size(3);
max_x = vol_size(4);
min_y = vol_size(1);
max_y = vol_size(2);

% Volume geometry object
vol_geom = astra_create_vol_geom(row_count, col_count, min_x, max_x, min_y, max_y);

% Normalize sino instead of volume
%% CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
sino = pi / 2 / length(angles) * sino;

%% Sinogram object
sino_id = astra_mex_data2d('create', '-sino', proj_geom, sino);

%% Volume object
vol_id = astra_mex_data2d('create', '-vol', vol_geom);

%% ASTRA config struct
switch lower( algorithm )
    case 'fbp'
        cfg = astra_struct('BP_CUDA');
    case 'fbp-astra'
        cfg = astra_struct('FBP_CUDA');
    case 'sirt'
        cfg = astra_struct('SIRT_CUDA');
    case 'sart'
        cfg = astra_struct('SART_CUDA');
    case 'cgls'
        cfg = astra_struct('CGLS_CUDA');
    case 'em'
        cfg = astra_struct('EM_CUDA');
end
if sum( strcmpi( algorithm, { 'sirt','sart' }) )
    if ~isempty( MinConstraint )
        cfg.option.MinConstraint = MinConstraint;
        
    end
    if ~isempty( MaxConstraint )
        cfg.option.MaxConstraint = MaxConstraint;
    end
end
cfg.ProjectionDataId = sino_id;
cfg.ReconstructionDataId = vol_id;

%% ASTRA create algorithm object from configuration struct
bp_id = astra_mex_algorithm('create', cfg);

%% Backprojection: iterate algorithm 
% SIRT_CUDA supports astra_mex_algorithm(‘get_res_norm’) to get the 2-norm
% of the difference between the projection data and the projection of the
% reconstruction.
if strcmpi( algorithm, 'fbp' )
    astra_mex_algorithm('iterate', bp_id, 1);
else
    astra_mex_algorithm('iterate', bp_id, iterations);
end
astra_mex_algorithm('delete', bp_id);
astra_mex_data2d('delete', sino_id)

%% Fetch data from ASTRA memory
vol = astra_mex_data2d('get_single', vol_id);
astra_mex_data2d('delete', vol_id)
%astra_mex_projector2d('clear')

% Required for adjoint?
%vol = pi / 2 / numel(angles) * vol;

end
