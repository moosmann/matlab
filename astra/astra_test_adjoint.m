function astra_test_adjoint()
% <A x, y> = <x, Ad y>

%% Paramters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% volume sizes
vol_shape_2d = [100, 100];
vol_shape_3d = [100, 100, 100];

% distances for fanflat and cone
source_origin = 300;
origin_detector = 100;

% parameters for 2D and 3D
det_col_count = 100; % for 3D use square detector
num_angles = 1;
det_width = 2 ;% for 3D use quadratic detector
angles = linspace2(0, 3*pi, num_angles); % excludes upper limit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mag = (source_origin+origin_detector)/source_origin;

% unit volume, unit projections
x = ones(vol_shape_2d);
y = ones(num_angles, det_col_count);

% Volume geometry
vol_geom = astra_create_vol_geom([vol_shape_2d(1), vol_shape_2d(2)]);

%% 2D parallel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% projection geometry
proj_geom = astra_create_proj_geom('parallel', det_width, det_col_count, angles);

%% line
d2_par_line.geom = proj_geom.type;
weight = 'line';
d2_par_line.weight = weight;
% FP
d2_par_line.Ax = astra_2d_fp(weight, proj_geom, vol_geom, x);
% BP
d2_par_line.Ady = astra_2d_bp(weight, proj_geom, vol_geom, y);
% inner products
d2_par_line.innProdProj = sum(d2_par_line.Ax(:) .* y(:));
d2_par_line.innProdVol = sum(x(:) .* d2_par_line.Ady(:));


%% linear
d2_par_linear.geom = proj_geom.type;
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
d2_par_strip.geom = proj_geom.type;
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
d2_par_cuda.geom = proj_geom.type;
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

% projection geometry
proj_geom = astra_create_proj_geom('fanflat', det_width, det_col_count, ...
    angles, source_origin, origin_detector);

%% linear
d2_fan_line.geom = proj_geom.type;
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
d2_fan_strip.geom = proj_geom.type;
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
d2_fan_cuda.geom = proj_geom.type;
weight = 'cuda';
d2_fan_cuda.weight = weight;
% FP
d2_fan_cuda.Ax = astra_2d_fp_cuda(proj_geom, vol_geom, x);
% BP
d2_fan_cuda.Ady = astra_2d_bp_cuda(proj_geom, vol_geom, y);
% inner products
d2_fan_cuda.innProdProj = sum(d2_fan_cuda.Ax(:) .* y(:)) * det_width/mag;
d2_fan_cuda.innProdVol = sum(x(:) .* d2_fan_cuda.Ady(:));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vol_shape_geom = [vol_shape_3d(2), vol_shape_3d(1), vol_shape_3d(3)];
det_row_count = det_col_count;

% unit volume, unit projections
x = ones(vol_shape_3d);
y = ones(det_row_count, num_angles, det_col_count);

% Volume geometry
vol_geom = astra_create_vol_geom(vol_shape_geom );

%% 3D parallel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% projection geometry
det_spacing_x = det_width;
det_spacing_y = det_width;
proj_geom = astra_create_proj_geom('parallel3d', det_spacing_x, ...
    det_spacing_y, det_row_count, det_col_count, angles);

%% cuda
d3_par_cuda.geom = proj_geom.type;
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

% projection geometry
det_spacing_x = det_width;
det_spacing_y = det_width;
proj_geom = astra_create_proj_geom('cone', det_spacing_x, det_spacing_y, ...
    det_row_count, det_col_count, angles, source_origin, origin_detector);

%% cuda
d3_cone_cuda.geom = proj_geom.type;
weight = 'cuda';
d3_cone_cuda.weight = weight;
% FP
d3_cone_cuda.Ax = astra_3d_fp_cuda(proj_geom, vol_geom, x);
% BP
d3_cone_cuda.Ady = astra_3d_bp_cuda(proj_geom, vol_geom, y);
% inner products
d3_cone_cuda.innProdProj = sum(d3_cone_cuda.Ax(:) .* y(:))  * det_width^2 / mag^2;
d3_cone_cuda.innProdVol = sum(x(:) .* d3_cone_cuda.Ady(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Print %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nnumber angles = %u, detector pixels width = %g', num_angles, det_width)
fprintf('\n2D: number voxels = (%u, %u), number detector pixels = %u', vol_shape_2d, det_col_count)
fprintf('\n3D: number voxels = (%u, %u, %u), detector pixels = (%u, %u)', vol_shape_3d, det_col_count, det_row_count)
fprintf('\nfanflat & cone: distance source orgin = %g, distance origin detector = %g, magnification = %g', source_origin, origin_detector, mag)

fprintf('\n')

pform = '\n%-24s';
pform2 = '%14g';
pform2s = '%14s';

p = @(format_specifier, f) fprintf(format_specifier, ... 
   f(d2_par_line), f(d2_par_linear), f(d2_par_strip), f(d2_par_cuda), f(d2_fan_line), f(d2_fan_strip), f(d2_fan_cuda), f(d3_par_cuda), f(d3_cone_cuda));

% geometry
fprintf(pform, 'geometry')
f = @(s) s.geom;
p(pform2s, f)
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
% mean
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

% inner product projection space
fprintf(pform, '<A x, y>/#pixel')
f = @(s) s.innProdProj/numel(s.Ax);
p(pform2, f)
% inner product volume space
fprintf(pform, '<x, Ad y>/#voxel')
f = @(s) s.innProdVol/numel(s.Ady);
p(pform2, f)


% ratio
fprintf(pform, '|<A x, y>/<x, Ad y>-1|')
f = @(s) abs(s.innProdProj / s.innProdVol-1);
p(pform2, f)
 % ratio
fprintf(pform, '|<x, Ad y>/<A x, y>-1|')
f = @(s) abs(s.innProdVol / s.innProdProj -1);
p(pform2, f)
% ratio
fprintf(pform, '<A x, y>/<x, Ad y>')
f = @(s) s.innProdProj / s.innProdVol;
p(pform2, f)
% ratio 
fprintf(pform, '<x, Ad y>/<A x, y>')
f = @(s) s.innProdVol / s.innProdProj;
p(pform2, f)

fprintf('\n')

%% Show %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
rows = 3;
cols = 3;

figure('Name','FP')
counter = 1;
for s = [d2_par_line, d2_par_linear, d2_par_strip, d2_par_cuda, d2_fan_line, d2_fan_strip, d2_fan_cuda, d3_par_cuda, d3_cone_cuda]
    subplot(rows ,cols,counter)
    if ndims(s.Ax) == 2
        imshow(s.Ax, [], 'InitialMagnification','fit')
    elseif ndims(s.Ax) == 3
        n = ceil(size(s.Ax, 3)/2);
        imshow(squeeze(s.Ax(:,:,n)), [], 'InitialMagnification','fit')
    end
    title(sprintf('%s %s', s.geom, s.weight), 'Interpreter', 'none')
    colorbar
    counter = counter + 1;
end

figure('Name','BP')
counter = 1;
for s = [d2_par_line, d2_par_linear, d2_par_strip, d2_par_cuda, d2_fan_line, d2_fan_strip, d2_fan_cuda, d3_par_cuda, d3_cone_cuda]
    subplot(rows ,cols,counter)
    if ndims(s.Ax) == 2
        imshow(s.Ady, [], 'InitialMagnification','fit')        
    elseif ndims(s.Ax) == 3
        n = ceil(size(s.Ady, 3)/2);
        imshow(squeeze(s.Ady(:,:,n)'), [], 'InitialMagnification','fit')
    end
    title(sprintf('%s %s', s.geom, s.weight), 'Interpreter', 'none')
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

% number angles = 5, detector pixels width = 2
% 2D: number voxels = (101, 101), number detector pixels = 100
% 3D: number voxels = (300, 300, 300), detector pixels = (100, 100)
% fanflat & cone: distance source orgin = 300, distance origin detector = 100, magnification = 1.33333
% 
% geometry                      parallel      parallel      parallel      parallel       fanflat       fanflat       fanflat    parallel3d          cone
%                                   line        linear         strip          cuda  line_fanflat strip_fanflat          cuda          cuda          cuda
% min(Ax)                              0             0             0             0             0             0             0       232.429       282.569
% max(Ax)                        124.843       124.843       249.687       249.686       127.436       194.979       127.436        370.82       400.917
% mean(Ax)                        50.906       50.9052        102.01        101.81       68.9897        102.01       68.9928       315.638        328.26
% min(Ad y)                            0             0       4.99998             5             0       4.99476             5             0             0
% max(Ad y)                      5.26344       4.41354       5.00002             5       5.83996             5             5             5             5
% mean(Ad y)                     2.49515       2.49511             5             5       3.38152       4.99999             5       2.33802       1.42656
% <A x, y>                         25453       25452.6         51005       50905.2       34494.9       51004.9       51744.6   6.31275e+07   3.69293e+07
% <x, Ad y>                        25453       25452.6         51005         51005       34494.9       51004.9         51005   6.31264e+07    3.8517e+07
% |<A x, y>/<x, Ad y>-1|     8.47551e-08   2.49664e-09   6.89464e-08    0.00195639    5.6175e-08   3.64336e-08     0.0145009   1.75865e-05     0.0412221
% |<x, Ad y>/<A x, y>-1|     8.47551e-08   2.49664e-09   6.89464e-08    0.00196023    5.6175e-08   3.64336e-08     0.0142937   1.75862e-05     0.0429945
% <A x, y>/<x, Ad y>                   1             1             1      0.998044             1             1        1.0145       1.00002      0.958778
% <x, Ad y>/<A x, y>                   1             1             1       1.00196             1             1      0.985706      0.999982       1.04299
