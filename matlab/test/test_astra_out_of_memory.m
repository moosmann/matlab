clear all
astra_clear
link_data = 1;
nxy = 2560;
nz = 1195;
num_proj = 5000;
angles = linspace2(0, pi, num_proj);
rao = 1178;
angles = double( angles );
algorithm = 'fbp';
astra_mex('set_gpu_index', 0:gpuDeviceCount() - 1);


sino = ones([nxy,num_proj,nz],'single');
vol_shape = [size( sino, 1), size( sino, 1), size(sino, 3) ];
vol_size = [];


%% Detector geometry
det_col_count = size( sino, 1);
det_row_count = size( sino, 3);
num_proj = size( sino, 2);


% Create geometry vector
vectors = zeros( numel(angles), 12);
tilt_camera = 0;
tilt_lamino = 0;
DetectorSpacingX = 1;
DetectorSpacingY = 1;
for nn = 1:num_proj
    theta = angles( nn );
    % source / ray direction
    vectors(nn,1) = + sin( theta );
    vectors(nn,2) = - cos( theta );
    vectors(nn,3) = - sin( tilt_lamino );

    % center of detector
    vectors(nn,4) = rao * cos( theta );
    vectors(nn,5) = rao * sin( theta );
    vectors(nn,6) = sin( tilt_lamino );

    % vector from detector pixel (0,0) to (0,1)
    vectors(nn,7) = cos( tilt_camera ) * cos( theta ) * DetectorSpacingX;
    vectors(nn,8) = cos( tilt_camera ) * sin( theta ) * DetectorSpacingX;
    vectors(nn,9) = cos( tilt_lamino ) * sin( tilt_camera ) * DetectorSpacingX;

    % vector from detector pixel (0,0) to (1,0)
    vectors(nn,10) = -sin( tilt_camera ) * cos( theta ) * DetectorSpacingY;
    vectors(nn,11) = -sin( tilt_camera ) * sin( theta ) * DetectorSpacingY;
    vectors(nn,12) = cos( tilt_lamino) * cos(tilt_camera) * DetectorSpacingY;
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
    case 'cgls'
        cfg = astra_struct('CGLS3D_CUDA');
end
if sum( strcmpi( algorithm, { 'sirt' }) )
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
tic;
if strcmpi( algorithm, 'fbp' )
    astra_mex_algorithm('iterate', bp_id, 1);
else
    astra_mex_algorithm('iterate', bp_id, iterations);
end
t = toc;
fprintf('\n duration: %f.0 s = %.1f min',t,t/60)
astra_mex_algorithm('delete', bp_id);
astra_mex_data3d('delete', sino_id)

%% Fetch data from ASTRA memory
if ~link_data
    vol = astra_mex_data3d('get_single', vol_id);
end
astra_mex_data3d('delete', vol_id)
astra_mex_projector3d('clear')


fprintf('\n')