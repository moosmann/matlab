% ASTRA toolbox reconstruction forwardard and backprojection for 3D vector
% geometry of a single slice. Test if angles project orthogonally.

clear all
aclear

crop = 2;

%% Parameter
N = 15;
Z = N - 1;
M = N + 1;
vol_shape_3d = [N, M, Z];
vol_shape_geom = [vol_shape_3d(2), vol_shape_3d(1), vol_shape_3d(3)];

angles = pi/2*[0];

det_row_count = Z;
det_col_count = N ; % for 3D use square detector
det_width = 1 ;% for 3D use quadratic detector

%% Phantom and projections
P = zeros(N,M,Z);
x = ceil((N-1)/4+1):floor(3/4*(N-1)+1);
y = ceil((M-1)/4+1):floor(3/4*(M-1)+1);
z = ceil((Z-1)/4+1):floor(3/4*(Z-1)+1);
P(x,y,z) = 1;
P = P + 0.0;

%% 3D parallel FP

% Volume geometry
vol_geom = astra_create_vol_geom(vol_shape_geom );

% projection geometry
det_spacing_x = det_width;
det_spacing_y = det_width;
proj_geom = astra_create_proj_geom('parallel3d', det_spacing_x, ...
    det_spacing_y, det_row_count, det_col_count, angles);

% create volume
vol_id = astra_mex_data3d('create','-vol', vol_geom, P);

% create sino
sino_id = astra_mex_data3d('create','-sino', proj_geom, 0);
    
% create struct
cfg = astra_struct('FP3D_CUDA');
cfg.ProjectionDataId = sino_id;
cfg.VolumeDataId = vol_id;
cfg.option.GPUindex = 1; 

% forward projection
alg_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('iterate', alg_id);
sino = astra_mex_data3d('get',sino_id);

%% Reco %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rotation_axis_offset = 0 - crop/2;
vol_shape = vol_shape_3d;
vol_size = [-vol_shape(1)/2, vol_shape(1)/2, -vol_shape(2)/2, vol_shape(2)/2, -vol_shape(3)/2, vol_shape(3)/2];
pixel_size = [det_width, det_width];
link_data = 0;
vol_shape = vol_shape .* [1 1 1];

%% unfiltered
tomo.angles = angles;
tomo.rot_axis_offset = rotation_axis_offset;
tomo.vol_shape = vol_shape;
tomo.vol_size = vol_size;
tomo.pixel_size = pixel_size;
tomo.link_data = link_data;
%sinoc = sino(1:end-crop,:,:);
sinoc = sino(1+crop:end,:,:);
bp = astra_parallel3D(tomo, sinoc);

figure('Name', sprintf('crop: %u', crop))
n = 1; m = 3;
subplot(n,m,1)
z2 = floor(Z/2);
imshow(P(:,:,z2), [],'InitialMagnification','fit'), title('phantom: P')

subplot(n,m,2)
imshow(squeeze(sinoc(:,1,:)), []), title('sino')

subplot(n,m,3)
imshow(bp(:,:,z2), [], 'InitialMagnification','fit'), title('bp')

whos P sino bp