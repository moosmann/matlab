% Test if back-projector 'Ad' is adjoint of forward projector 'A':
% <A x, y> = <x, Ad y>
clear
astra_clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2D parallel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vol_shape = [100, 100];
vol_shape_geom = [vol_shape(2), vol_shape(1)];
det_col_count = 100;
num_angles = 200;
det_width = 1.3;
angles = linspace2(0, pi, num_angles); % excludes upper limit

fprintf('\n2D PARALLEL: #vox = (%u, %u), #ang = %u, #pix = %u, d_pix = %g', vol_shape, num_angles, det_col_count, det_width)
% unit volume, unit projections
x = ones(vol_shape);
y = ones(num_angles, det_col_count);

% Volume geometry
vol_geom = astra_create_vol_geom(vol_shape_geom);

% projection geometry
proj_geom = astra_create_proj_geom('parallel', det_width, det_col_count, angles);

%% FB
% create volume
vol_id = astra_mex_data2d('create','-vol', vol_geom, x);

% create sino
sino_id = astra_mex_data2d('create','-sino', proj_geom, 0);
    
% create struct
cfg = astra_struct('FP_CUDA');
cfg.ProjectionDataId = sino_id;
cfg.VolumeDataId = vol_id;
cfg.option.GPUindex = 1; % 'GeForce GTX 980'

% create forward projection
alg_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('iterate', alg_id);
Ax = astra_mex_data2d('get',sino_id);
Axpara = Ax;

astra_clear

%% BP
% create volume
vol_id = astra_mex_data2d('create','-vol', vol_geom, 0);

% store sino
sino_id = astra_mex_data2d('create','-sino', proj_geom, y);

% create struct
cfg = astra_struct('BP_CUDA');
cfg.ProjectionDataId = sino_id;
cfg.ReconstructionDataId = vol_id;
cfg.option.GPUindex = 1; % 'GeForce GTX 980'

% create back projection
alg_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('iterate', alg_id);
Ady = astra_mex_data2d('get', vol_id);

%% Check adjoint: <A x, y> = <x, Ad y>

l = sum(Ax(:) .* y(:));
r = sum(x(:) .* Ady(:));

fprintf('\n %-20s | %-20s | %-20s | %-20s', '<A x, y>', '<x, Ad y>', '<A x, y>/<x, Ad y> ', '<x, Ad y>/<A x, y>')
fprintf('\n %20g | %20g | %20.8g | %20.8g', l, r, l/r, r/l)

subplot(4,2,1)
imshow(Ax, [], 'InitialMagnification','fit')
title('2D parallel FP')
subplot(4,2,2)
imshow(Ady, [], 'InitialMagnification','fit')
title('2D parallel BP')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2D fanflat %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

astra_clear
%clear

vol_shape = [101, 103];
det_col_count = 99;
num_angles = 111;
det_width = 0.3;
source_origin = 100;
origin_detector = 30;
angles = linspace2(0, 2*pi, num_angles); % excludes upper limit
mag = (source_origin+origin_detector)*vol_shape(1)/source_origin/vol_shape(1);
l_fac = det_width;

fprintf('\n2D FANFLAT :')
fprintf('#vox = (%u, %u), #ang = %u, #pix = %u, d_pix = %g', vol_shape, num_angles, det_col_count, det_width)
fprintf(', d_so = %g, d_od = %g, mag = %g', source_origin, origin_detector, mag)

% unit volume, unit projections
x = ones(vol_shape);
y = ones(num_angles, det_col_count);
%y(1, :) = 1:det_col_count;

% Volume geometry
vol_geom = astra_create_vol_geom(vol_shape);

% projection geometry
proj_geom = astra_create_proj_geom('fanflat', det_width, det_col_count, ...
    angles, source_origin, origin_detector);

%% FP
% create volume
vol_id = astra_mex_data2d('create','-vol', vol_geom, x);

% create sino
sino_id = astra_mex_data2d('create','-sino', proj_geom, 0);

% create struct
cfg = astra_struct('FP_CUDA');
cfg.ProjectionDataId = sino_id;
cfg.VolumeDataId = vol_id;
cfg.option.GPUindex = 1; % 'GeForce GTX 980'

% create forward projection
alg_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('iterate', alg_id);
Ax = astra_mex_data2d('get', sino_id);

%% BP
% store sino
sino_id = astra_mex_data2d('create','-sino', proj_geom, y);

% create struct
cfg = astra_struct('BP_CUDA');
cfg.ProjectionDataId = sino_id;
cfg.ReconstructionDataId = vol_id;
cfg.option.GPUindex = 1; % 'GeForce GTX 980'

% create sinogram
alg_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('iterate', alg_id);
Ady = astra_mex_data2d('get', vol_id);

%% Check adjoint: <A x, y> = <x, Ad y>

l = sum(Ax(:) .* y(:)) * l_fac;
r = sum(x(:) .* Ady(:));

fprintf('\n %-20s | %-20s | %-20s | %-20s', '<A x, y>', '<x, Ad y>', '<A x, y>/<x, Ad y> ', '<x, Ad y>/<A x, y>')
fprintf('\n %20g | %20g | %20.8g | %20.8g', l, r, l/r, r/l)

subplot(4,2,3)
imshow(Ax, [], 'InitialMagnification','fit')
title('2D fanflat FP')
subplot(4,2,4)
imshow(Ady, [], 'InitialMagnification','fit')
title('2D fanflat BP')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3D parallel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vol_shape = [101, 101, 101];
vol_shape_geom = [vol_shape(2), vol_shape(1), vol_shape(3)];
det_col_count = 100;
det_row_count = det_col_count;
num_angles = 200;
det_width = 1.3;
angles = linspace2(0, pi, num_angles); % excludes upper limit

fprintf('\n3D PARALLEL: ')
fprintf('#vox = (%u, %u, %u), #ang = %u', vol_shape, num_angles)
fprintf(', #pix = (%u, %u), d_pix = %g', det_col_count, det_row_count, det_width)

% unit volume, unit projections
x = ones(vol_shape);
y = ones(det_row_count, num_angles, det_col_count);

% Volume geometry
vol_geom = astra_create_vol_geom(vol_shape_geom );

% projection geometry
det_spacing_x = det_width;
det_spacing_y = det_width;
proj_geom = astra_create_proj_geom('parallel3d', det_spacing_x, ...
    det_spacing_y, det_row_count, det_col_count, angles);

%% FB
% create volume
vol_id = astra_mex_data3d('create','-vol', vol_geom, x);

% create sino
sino_id = astra_mex_data3d('create','-sino', proj_geom, 0);
    
% create struct
cfg = astra_struct('FP3D_CUDA');
cfg.ProjectionDataId = sino_id;
cfg.VolumeDataId = vol_id;
cfg.option.GPUindex = 1; % 'GeForce GTX 980'

% forward projection
alg_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('iterate', alg_id);
Ax = astra_mex_data3d('get',sino_id);
Axpara3 = Ax;

astra_clear

%% BP
% create volume
vol_id = astra_mex_data3d('create','-vol', vol_geom, 0);

% store sino
sino_id = astra_mex_data3d('create','-sino', proj_geom, y);

% create struct
cfg = astra_struct('BP3D_CUDA');
cfg.ProjectionDataId = sino_id;
cfg.ReconstructionDataId = vol_id;
cfg.option.GPUindex = 1; % 'GeForce GTX 980'

% create back projection
alg_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('iterate', alg_id);
Ady = astra_mex_data3d('get', vol_id);

astra_clear

%% Check adjoint: <A x, y> = <x, Ad y>

l = sum(Ax(:) .* y(:)) * det_width^2;
%l = sum(Ax(:) .* y(:));
r = sum(x(:) .* Ady(:));

fprintf('\n %-20s | %-20s | %-20s | %-20s', '<A x, y>', '<x, Ad y>', '<A x, y>/<x, Ad y> ', '<x, Ad y>/<A x, y>')
fprintf('\n %20g | %20g | %20.8g | %20.8g', l, r, l/r, r/l)

subplot(4,2,5)
imshow(Ax(:,:,ceil(size(Ax,3)/2)), [], 'InitialMagnification','fit')
%imshow(squeeze(Ax(ceil(size(Ax,1)/2),:,:)), [], 'InitialMagnification','fit')
title('3D parallel FP')
subplot(4,2,6)
imshow(Ady(:,:,ceil(size(Ax,3)/2)), [], 'InitialMagnification','fit')
%imshow(squeeze(Ady(ceil(size(Ax,1)/2),:,:)), [], 'InitialMagnification','fit')
title('3D parallel BP')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3D cone %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

astra_clear
% clear

vol_shape = [101, 101, 101];
vol_shape_geom = [vol_shape(2), vol_shape(1), vol_shape(3)];
det_col_count = 101;
det_row_count = 101;
num_angles = 1;
det_width = 1;
source_origin = 51;
origin_detector = 10000;
angles = linspace2(0, 2*pi, num_angles); % excludes upper limit
mag = (source_origin+origin_detector)*vol_shape(1)/source_origin/vol_shape(1);

fprintf('\n3D CONE:')
fprintf(' #vox = (%u, %u, %u), #ang = %u, #pix = (%u, %u), d_pix = %g', vol_shape, num_angles, det_col_count, det_row_count, det_width)
fprintf(', d_so = %g, d_od = %g, mag = %g', source_origin, origin_detector, mag)

% unit volume, unit projections
x = ones(vol_shape);
y = ones(det_col_count, num_angles, det_row_count);

% Volume geometry
vol_geom = astra_create_vol_geom(vol_shape_geom);

% projection geometry
det_spacing_x = det_width;
det_spacing_y = det_width;
proj_geom = astra_create_proj_geom('cone', det_spacing_x, det_spacing_y, ...
    det_row_count, det_col_count, angles, source_origin, origin_detector);

astra_clear

%% FP
% create volume
vol_id = astra_mex_data3d('create','-vol', vol_geom, x);

% create sino
sino_id = astra_mex_data3d('create','-sino', proj_geom, 0);
    
% create struct
cfg = astra_struct('FP3D_CUDA');
cfg.ProjectionDataId = sino_id;
cfg.VolumeDataId = vol_id;
cfg.option.GPUindex = 1; % 'GeForce GTX 980'

% create forward projection
alg_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('iterate', alg_id);
Ax = astra_mex_data3d('get',sino_id);

%% BP of y
% create volume
vol_id = astra_mex_data3d('create','-vol', vol_geom, 0);

% store sino
sino_id = astra_mex_data3d('create','-sino', proj_geom, y);
    
% create struct
cfg = astra_struct('BP3D_CUDA');
cfg.ProjectionDataId = sino_id;
cfg.ReconstructionDataId = vol_id;
cfg.option.GPUindex = 1; % 'GeForce GTX 980'

% create backprojection
alg_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('iterate', alg_id);
Ady = astra_mex_data3d('get',vol_id);

astra_clear

%% BP of Ax
% create volume
vol_id = astra_mex_data3d('create','-vol', vol_geom, 0);

% store sino
sino_id = astra_mex_data3d('create','-sino', proj_geom, Ax);
    
% create struct
cfg = astra_struct('BP3D_CUDA');
cfg.ProjectionDataId = sino_id;
cfg.ReconstructionDataId = vol_id;
cfg.option.GPUindex = 1; % 'GeForce GTX 980'

% create backprojection
alg_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('iterate', alg_id);
AdAx = astra_mex_data3d('get',vol_id);


astra_clear

%% Check adjoint: <A x, y> = <x, Ad y>

l = sum(Ax(:) .* y(:)) * det_width^2;
r = sum(x(:) .* Ady(:));
l1 = sum(Ax(:).^2);
r1 = sum(x(:) .* AdAx(:));

fprintf('\n %-15s | %-15s | %-15s | %-15s | %-15s | %-15s', '<A x, y>', '<x, Ad y>', '<A x, y>/<x, Ad y> ', '<x, Ad y>/<A x, y>')
fprintf('\n %15g | %15g | %15.8g | %15.8g | %-15.8g | %-15.8g', l, r, l/r, r/l, l1/r1, r1/l1)

subplot(4,2,7)
imshow(Ax(:,:,ceil(size(Ax,3)/2)), [], 'InitialMagnification','fit')
%imshow(squeeze(Ax(ceil(size(Ax,1)/2),:,:)), [], 'InitialMagnification','fit')
title('3D cone FP')
subplot(4,2,8)
imshow(Ady(:,:,ceil(size(Ax,3)/2)), [], 'InitialMagnification','fit')
%imshow(squeeze(Ady(ceil(size(Ax,1)/2),:,:)), [], 'InitialMagnification','fit')
title('3D cone BP')

fprintf('\n')