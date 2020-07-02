% ASTRA toolbox reconstruction forwardard and backprojection for 3D vector
% geometry of a single slice. Test if angles project orthogonally.

clear all
aclear

%% Parameter
N = 8;
M = N +0;
vol_shape_3d = [N, M, 1];
vol_shape_geom = [vol_shape_3d(2), vol_shape_3d(1), vol_shape_3d(3)];

angles = pi/2*[0 1 2 3 ];

det_row_count = 1;
det_col_count = N + 100; % for 3D use square detector
det_width = 1 ;% for 3D use quadratic detector

%% Phantom and projections
P = zeros(N,M);
x = ceil((N-1)/4+1):floor(3/4*(N-1)+1);
y = ceil((M-1)/4+1):floor(3/4*(M-1)+1);
P(x,y) = 1;
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
rotation_axis_offset = 0;
vol_shape = vol_shape_3d;
vol_size = [-vol_shape(1)/2, vol_shape(1)/2, -vol_shape(2)/2, vol_shape(2)/2, -vol_shape(3)/2, vol_shape(3)/2];
pixel_size = [det_width, det_width];
link_data = 0;
vol_shape = vol_shape .* [1 1 1];

%% unfiltered
bp = astra_parallel3D(sino, angles, rotation_axis_offset, vol_shape, vol_size, pixel_size, link_data);

figure(1)
%subplot(2,2,1)
imshow(P, [],'InitialMagnification','fit'), title('phantom: P')
%subplot(2,2,2)
figure(2)
imshow(bp, [], 'InitialMagnification','fit'), title('bp')
%subplot(2,2,3:4), imshow(sino', []), title('sino')
