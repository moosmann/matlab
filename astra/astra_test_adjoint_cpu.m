function astra_test_adjoint_cpu()
% <A x, y> = <x, Ad y>

%% Paramters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% volume sizes
vol_shape_2d = [101, 101];
vol_shape_3d = [101, 101, 101];

% distances for fanflat and cone
source_origin = 100000;
origin_detector = 100;

% parameters for 2D and 3D
det_col_count = 99; % for 3D use square detector
num_angles = 1;
det_width = 0.3;
angles = linspace2(0, pi, num_angles); % excludes upper limit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vol_shape = vol_shape_2d;

% unit volume, unit projections
x = ones(vol_shape);
y = ones(num_angles, det_col_count);

%% 2D parallel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n2D parallel: #vox = (%u, %u), #ang = %u, #pix = %u, d_pix = %g', vol_shape, num_angles, det_col_count, det_width)

% Volume geometry
vol_geom = astra_create_vol_geom(vol_shape);

% projection geometry
proj_geom = astra_create_proj_geom('parallel', det_width, det_col_count, angles);

%% linear
weight = 'linear'; % Josephson kernel
d2_par_linear.weight = weight;
% FP
d2_par_linear.Ax = astra_2d_fp(weight, proj_geom, vol_geom, x);
% BP
d2_par_linear.Ady = astra_2d_bp(weight, proj_geom, vol_geom, y);
% inner products
d2_par_linear.innProdProj = sum(d2_par_linear.Ax(:) .* y(:));
d2_par_linear.innProdVol = sum(x(:) .* d2_par_linear.Ady(:));

%% strip
weight = 'strip';
d2_par_strip.weight = weight;
% FP
d2_par_strip.Ax = astra_2d_fp(weight, proj_geom, vol_geom, x);
% BP
d2_par_strip.Ady = astra_2d_bp(weight, proj_geom, vol_geom, y);
% inner products
d2_par_strip.innProdProj = sum(d2_par_strip.Ax(:) .* y(:));
d2_par_strip.innProdVol = sum(x(:) .* d2_par_strip.Ady(:));

%% cuda
weight = 'cuda';
d2_par_cuda.weight = weight;
% FP
d2_par_cuda.Ax = astra_2d_fp_cuda(proj_geom, vol_geom, x);
% BP
d2_par_cuda.Ady = astra_2d_bp_cuda(proj_geom, vol_geom, y);
% inner products
d2_par_cuda.innProdProj = sum(d2_par_cuda.Ax(:) .* y(:));
d2_par_cuda.innProdVol = sum(x(:) .* d2_par_cuda.Ady(:));


%% 2D fanflat %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mag2d = (source_origin+origin_detector)*vol_shape(1)/source_origin/vol_shape(1);

fprintf('\n2D fanflat :')
fprintf(' #vox = (%u, %u), #ang = %u, #pix = %u, d_pix = %g', vol_shape, num_angles, det_col_count, det_width)
fprintf(', d_so = %g, d_od = %g, mag = %g', source_origin, origin_detector, mag2d)

% Volume geometry
vol_geom = astra_create_vol_geom(vol_shape);

% projection geometry
proj_geom = astra_create_proj_geom('fanflat', det_width, det_col_count, ...
    angles, source_origin, origin_detector);

%% linear
weight = 'line_fanflat';
d2_fan_line.weight = weight;
% FP
d2_fan_line.Ax = astra_2d_fp(weight, proj_geom, vol_geom, x);
% BP
d2_fan_line.Ady = astra_2d_bp(weight, proj_geom, vol_geom, y);
% inner products
d2_fan_line.innProdProj = sum(d2_fan_line.Ax(:) .* y(:));
d2_fan_line.innProdVol = sum(x(:) .* d2_fan_line.Ady(:));

%% strip
weight = 'strip_fanflat';
d2_fan_strip.weight = weight;
% FP
d2_fan_strip.Ax = astra_2d_fp(weight, proj_geom, vol_geom, x);
% BP
d2_fan_strip.Ady = astra_2d_bp(weight, proj_geom, vol_geom, y);
% inner products
d2_fan_strip.innProdProj = sum(d2_fan_strip.Ax(:) .* y(:));
d2_fan_strip.innProdVol = sum(x(:) .* d2_fan_strip.Ady(:));

%% cuda
weight = 'cuda';
d2_fan_cuda.weight = weight;
% FP
d2_fan_cuda.Ax = astra_2d_fp_cuda(proj_geom, vol_geom, x);
% BP
d2_fan_cuda.Ady = astra_2d_bp_cuda(proj_geom, vol_geom, y);
% inner products
d2_fan_cuda.innProdProj = sum(d2_fan_cuda.Ax(:) .* y(:)) * det_width;
d2_fan_cuda.innProdVol = sum(x(:) .* d2_fan_cuda.Ady(:));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vol_shape = vol_shape_3d;
vol_shape_geom = [vol_shape(2), vol_shape(1), vol_shape(3)];
det_row_count = det_col_count;

% unit volume, unit projections
x = ones(vol_shape);
y = ones(det_row_count, num_angles, det_col_count);

%% 3D parallel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n3D parallel:')
fprintf(' #vox = (%u, %u, %u), #ang = %u', vol_shape, num_angles)
fprintf(', #pix = (%u, %u), d_pix = %g', det_col_count, det_row_count, det_width)

% Volume geometry
vol_geom = astra_create_vol_geom(vol_shape_geom );

% projection geometry
det_spacing_x = det_width;
det_spacing_y = det_width;
proj_geom = astra_create_proj_geom('parallel3d', det_spacing_x, ...
    det_spacing_y, det_row_count, det_col_count, angles);

%% cuda
weight = 'cuda';
d3_par_cuda.weight = weight;
% FP
d3_par_cuda.Ax = astra_3d_fp_cuda(proj_geom, vol_geom, x);
% BP
d3_par_cuda.Ady = astra_3d_bp_cuda(proj_geom, vol_geom, y);
% inner products
d3_par_cuda.innProdProj = sum(d3_par_cuda.Ax(:) .* y(:)) * det_width^2;
d3_par_cuda.innProdVol = sum(x(:) .* d3_par_cuda.Ady(:));

%% cone %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mag3d = (source_origin+origin_detector)*vol_shape(1)/source_origin/vol_shape(1);

fprintf('\n3D cone:')
fprintf('     #vox = (%u, %u, %u), #ang = %u, #pix = (%u, %u), d_pix = %g', vol_shape, num_angles, det_col_count, det_row_count, det_width)
fprintf(', d_so = %g, d_od = %g, mag = %g', source_origin, origin_detector, mag3d)

% Volume geometry
vol_geom = astra_create_vol_geom(vol_shape_geom );

% projection geometry
det_spacing_x = det_width;
det_spacing_y = det_width;
proj_geom = astra_create_proj_geom('cone', det_spacing_x, det_spacing_y, ...
    det_row_count, det_col_count, angles, source_origin, origin_detector);

%% cuda
weight = 'cuda';
d3_cone_cuda.weight = weight;
% FP
d3_cone_cuda.Ax = astra_3d_fp_cuda(proj_geom, vol_geom, x);
% BP
d3_cone_cuda.Ady = astra_3d_bp_cuda(proj_geom, vol_geom, y);
% inner products
d3_cone_cuda.innProdProj = sum(d3_cone_cuda.Ax(:) .* y(:))  * det_width^2;
d3_cone_cuda.innProdVol = sum(x(:) .* d3_cone_cuda.Ady(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Print %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n')

pform = '\n%-24s';
pform2 = '%14g';
pform2s = '%14s';
pform2sr = '%-14s';

fprintf(pform, '')
fprintf(pform2sr, ' 2D parallel', '', '', ' 2D fanflat', '', '', ' 3D parallel', ' 3D cone')

p = @(format_specifier, f) fprintf(format_specifier, ... 
   f(d2_par_linear), f(d2_par_strip), f(d2_par_cuda), f(d2_fan_line), f(d2_fan_strip), f(d2_fan_cuda), f(d3_par_cuda), f(d3_cone_cuda));

% weights
fprintf(pform, '')
f = @(s) s.weight;
p(pform2s, f)

% min
fprintf(pform, 'min(Ax)')
f = @(s) min(s.Ax(:));
p(pform2, f)

% max
fprintf(pform, 'max(Ax)')
f = @(s) max(s.Ax(:));
p(pform2, f)

% max
fprintf(pform, 'mean(Ax)')
f = @(s) mean(s.Ax(:));
p(pform2, f)


% min
fprintf(pform, 'min(Ad y)')
f = @(s) min(s.Ady(:));
p(pform2, f)

% max
fprintf(pform, 'max(Ad y)')
f = @(s) max(s.Ady(:));
p(pform2, f)

% max
fprintf(pform, 'mean(Ad y)')
f = @(s) mean(s.Ady(:));
p(pform2, f)

% inner product projection space
fprintf(pform, '<A x, y>')
f = @(s) s.innProdProj;
p(pform2, f)

% inner product volume space
fprintf(pform, '<x, Ad y>')
f = @(s) s.innProdVol;
p(pform2, f)

% ratio
fprintf(pform, '|<A x, y>/<x, Ad y>-1|')
f = @(s) abs(s.innProdProj / s.innProdVol-1);
p(pform2, f)
 
fprintf(pform, '|<x, Ad y>/<A x, y>-1|')
f = @(s) abs(s.innProdVol / s.innProdProj -1);
p(pform2, f)

fprintf(pform, '<A x, y>/<x, Ad y>')
f = @(s) s.innProdProj / s.innProdVol;
p(pform2, f)
 
fprintf(pform, '<x, Ad y>/<A x, y>')
f = @(s) s.innProdVol / s.innProdProj;
p(pform2, f)

fprintf('\n')

%% Show %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
rows = 2;
cols = 3;

figure('Name','FP')
counter = 1;
for s = [d2_par_linear, d2_par_strip, d2_par_cuda, d2_fan_line, d2_fan_strip, d2_fan_cuda]
    subplot(rows ,cols,counter)    
    imshow(s.Ax, [], 'InitialMagnification','fit')
    title(sprintf('%s', s.weight), 'Interpreter', 'none')
    colorbar
    counter = counter + 1;
end

figure('Name','BP')
counter = 1;
for s = [d2_par_linear, d2_par_strip, d2_par_cuda, d2_fan_line, d2_fan_strip, d2_fan_cuda]
    subplot(rows ,cols,counter)    
    imshow(s.Ady, [], 'InitialMagnification','fit')
    title(sprintf('%s', s.weight), 'Interpreter', 'none')
    colorbar
    counter = counter + 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fp = astra_2d_fp(proj_weight, proj_geom, vol_geom, vol)

% projector
proj_id = astra_create_projector(proj_weight, proj_geom, vol_geom);

% create volume
vol_id = astra_mex_data2d('create','-vol', vol_geom, vol);

% create sino
sino_id = astra_mex_data2d('create','-sino', proj_geom, 0);
    
% create struct
cfg = astra_struct('FP');
cfg.ProjectorId = proj_id;
cfg.ProjectionDataId = sino_id;
cfg.VolumeDataId = vol_id;

% create forward projection
alg_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('iterate', alg_id);
fp = astra_mex_data2d('get',sino_id);

astra_clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bp = astra_2d_bp(proj_weight, proj_geom, vol_geom, proj)

% projector
proj_id = astra_create_projector(proj_weight, proj_geom, vol_geom);

% create volume
vol_id = astra_mex_data2d('create','-vol', vol_geom, 0);

% store sino
sino_id = astra_mex_data2d('create','-sino', proj_geom, proj);

% create struct
cfg = astra_struct('BP');
cfg.ProjectorId = proj_id;
cfg.ProjectionDataId = sino_id;
cfg.ReconstructionDataId = vol_id;

% create back projection
alg_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('iterate', alg_id);
bp = astra_mex_data2d('get', vol_id);

astra_clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fp = astra_2d_fp_cuda(proj_geom, vol_geom, vol)
% 2D CUDA forward projector. Caution: fanflat geomatry requires
% multiplication of projection by detector pixel width, which is not the
% case for parallel geometries or when using CUDA.

% create volume
vol_id = astra_mex_data2d('create','-vol', vol_geom, vol);

% create sino
sino_id = astra_mex_data2d('create','-sino', proj_geom, 0);
    
% create struct
cfg = astra_struct('FP_CUDA');
cfg.ProjectionDataId = sino_id;
cfg.VolumeDataId = vol_id;

% create forward projection
alg_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('iterate', alg_id);
fp = astra_mex_data2d('get',sino_id);

% if strcmp(proj_geom.type,'fanflat')
%     fp = proj_geom.DetectorWidth * fp;
% end

astra_clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bp = astra_2d_bp_cuda(proj_geom, vol_geom, proj)

% create volume
vol_id = astra_mex_data2d('create','-vol', vol_geom, 0);

% store sino
sino_id = astra_mex_data2d('create','-sino', proj_geom, proj);

% create struct
cfg = astra_struct('BP_CUDA');
cfg.ProjectionDataId = sino_id;
cfg.ReconstructionDataId = vol_id;

% create back projection
alg_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('iterate', alg_id);
bp = astra_mex_data2d('get', vol_id);

astra_clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fp = astra_3d_fp_cuda(proj_geom, vol_geom, vol)

% create volume
vol_id = astra_mex_data3d('create','-vol', vol_geom, vol);

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
fp = astra_mex_data3d('get',sino_id);

%fp = fp * proj_geom.DetectorSpacingX * proj_geom.DetectorSpacingY;

astra_clear

%% BP
function bp = astra_3d_bp_cuda(proj_geom, vol_geom, proj)

% create volume
vol_id = astra_mex_data3d('create','-vol', vol_geom, 0);

% store sino
sino_id = astra_mex_data3d('create','-sino', proj_geom, proj);

% create struct
cfg = astra_struct('BP3D_CUDA');
cfg.ProjectionDataId = sino_id;
cfg.ReconstructionDataId = vol_id;
cfg.option.GPUindex = 1; % 'GeForce GTX 980'

% create back projection
alg_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('iterate', alg_id);
bp = astra_mex_data3d('get', vol_id);

astra_clear
