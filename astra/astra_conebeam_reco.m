function [rec, sino] = astra_conebeam_reco(DataSetNum,NumIterations,PadHor_PadVer,scaleFactor,lambda)

if nargin < 1
    DataSetNum = 10;
end
if nargin < 2
    NumIterations = 6;
end
if nargin < 3
    PadHor_PadVer = [0 0];
end
if nargin < 4
    scaleFactor = 1;
end
if nargin < 5
    lambda = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
astra_clear
astra_set_gpu_index(0)

padding.size_hor = PadHor_PadVer(1);
padding.size_ver = PadHor_PadVer(2);
padding.valueOrMethod = 'replicate';

%% Data sets

% sino matrix: the object is initialized with the contents of this matrix.
% The matrix must be of size (u,angles,v), where u is the number of columns
% of the detector and v the number of rows as defined in the projection
% geometry 

data_dir = [ getenv('HOME') '/data/gate/'];

nn = 1;
data(nn).dir = data_dir ;
data(nn).filename = 'detector_two_spheres_Astra_20150313';
data(nn).fieldname = 'detector_astra';
data(nn).permuteOrder = [3 2 1];
data(nn).fullAngle_rad = - 2 * pi;
data(nn).source_origin_mm = 280; 
data(nn).origin_det_mm = 20;
data(nn).detector_width_mm = 50;

nn = 2;
data(nn).dir = data_dir ;
data(nn).filename = '20150317_water_spheres_all_photons';
data(nn).fieldname = 'detector_astra';
data(nn).permuteOrder = [1 2 3];
data(nn).fullAngle_rad = - 2 * pi;
data(nn).source_origin_mm = 280; 
data(nn).origin_det_mm = 20;
data(nn).detector_width_mm = 50;

nn = 3;
data(nn).dir = data_dir ;
data(nn).filename = '20150317_water_spheres_primary_photons';
data(nn).fieldname = 'detector_astra';
data(nn).permuteOrder = [1 2 3];
data(nn).fullAngle_rad = - 2 * pi;
data(nn).source_origin_mm = 280; 
data(nn).origin_det_mm = 20;
data(nn).detector_width_mm = 50;

nn = 4;
data(nn).dir = data_dir ;
data(nn).filename = '20150317_water_spheres_all_photons_high_act';
data(nn).fieldname = 'detector_astra';
data(nn).permuteOrder = [1 2 3];
data(nn).fullAngle_rad = - 2 * pi;
data(nn).source_origin_mm = 280; 
data(nn).origin_det_mm = 20;
data(nn).detector_width_mm = 50;

nn = 5;
data(nn).dir = data_dir ;
data(nn).filename = '20150320_central_water_spheres';
data(nn).fieldname = 'detector_astra';
data(nn).permuteOrder = [1 2 3];
data(nn).fullAngle_rad = - 2 * pi;
data(nn).source_origin_mm = 190; 
data(nn).origin_det_mm = 110;
data(nn).detector_width_mm = 50;

nn = 6;
data(nn).dir = data_dir ;
data(nn).filename = '20150320_shifted_water_spheres';
data(nn).fieldname = 'detector_astra';
data(nn).permuteOrder = [1 2 3];
data(nn).fullAngle_rad = - 2 * pi;
data(nn).source_origin_mm = 190; 
data(nn).origin_det_mm = 110;
data(nn).detector_width_mm = 50;

nn = 7;
data(nn).dir = data_dir ;
data(nn).filename = '20150324_water_spheres_larger_det';
data(nn).fieldname = 'detector_astra';
data(nn).permuteOrder = [1 2 3];
data(nn).fullAngle_rad = - 2 * pi;
data(nn).source_origin_mm = 190; 
data(nn).origin_det_mm = 110;
data(nn).detector_width_mm = 100;

nn = 8;
data(nn).dir = data_dir ;
data(nn).filename = '20150330_voxelbox_30x30';
data(nn).fieldname = 'detector_astra';
data(nn).permuteOrder = [1 2 3];
data(nn).fullAngle_rad = - 2 * pi;
data(nn).source_origin_mm = 190; 
data(nn).origin_det_mm = 110;
data(nn).detector_width_mm = 100;

nn = 9;
data(nn).dir = data_dir ;
data(nn).filename = '20150410_water_spheres_high_act';
data(nn).fieldname = 'detector_astra';
data(nn).permuteOrder = [1 2 3];
data(nn).fullAngle_rad =  - 2 * pi;
data(nn).source_origin_mm = 190; 
data(nn).origin_det_mm = 110;
data(nn).detector_width_mm = 100;

nn = 10;
data(nn).dir = [ getenv('HOME') '/data/matlab/'];
data(nn).filename = 'conebeam_2balls_h100_p360_v100';
data(nn).fieldname = 'proj_data';
data(nn).permuteOrder = [1 2 3];
data(nn).fullAngle_rad = 2 * pi;
data(nn).source_origin_mm = 190; 
data(nn).origin_det_mm = 110;
data(nn).detector_width_mm = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scale_fac = scaleFactor(1);

% Read data
data = data(DataSetNum);
sino = load( [data.dir data.filename '.mat'] );
%sino = single( getfield( sino, data.fieldname ) );
sino = single (sino.(data.fieldname) );
sino = permute( sino, data.permuteOrder );

num_proj = size(sino,2);
num_pixel_hor_0 = size(sino,1);
num_pixel_ver_0 = size(sino,3);

fprintf('\nData set: %s', data.filename )
fprintf('\nRaw projections: [Min, Max] = [%g, %g]', min( sino(:) ), max( sino(:) ) )

% Normalise
if DataSetNum ~= 10
    sino = sino + 1e1;
    sino =  1 - log( (sino) / max( sino(:)) ) ;
    sino = 1 * normat( sino );
end
fprintf('\nProjections: [Min, Max] = [%g, %g]', min( sino(:) ), max( sino(:) ) )

% Padding
if padding.size_hor < 0
    padding.size_hor = abs(padding.size_hor) * ceil( ( sqrt(2) -1 ) / 2 * num_pixel_hor_0 );    
end
if padding.size_ver < 0
    padding.size_ver = abs(padding.size_ver) * ceil( ( sqrt(2) -1 ) / 2 * num_pixel_ver_0 );    
end
sino = padarray(sino,[padding.size_hor 0 padding.size_ver],padding.valueOrMethod,'both');

%num_proj = size(sino,2);
num_pixel_hor  = size(sino,1);
num_pixel_ver  = size(sino,3);

fprintf( '\nDetector pixels [h x v]: %4u x %4u (original), %4u %4u (padded)', num_pixel_hor_0, num_pixel_ver_0, num_pixel_hor, num_pixel_ver )
fprintf( '\nNumber of projection: %4u', num_proj)


num_voxel_hor = scale_fac * num_pixel_hor_0 ;
num_voxel_ver = scale_fac * num_pixel_ver_0 ;
%x = round( 1 + (num_voxel_hor-1) * [x 1-x] );x = x(1):x(2);


%% Parameters in physical units
% length in mm 
det_width_mm = data.detector_width_mm ;
pixel_width_hor_mm = det_width_mm / num_pixel_hor_0 ;
pixel_width_ver_mm = det_width_mm / num_pixel_ver_0 ;
source_origin_mm = data.source_origin_mm;
origin_det_mm = data.origin_det_mm;
%source_det_mm = source_origin_mm + origin_det_mm;
%vol_width_mm = 2 * source_origin_mm * det_width_mm / ( 2 * source_det_mm + det_width_mm );
vol_width_mm = det_width_mm;
voxel_width_mm = vol_width_mm / num_voxel_hor;
fprintf( '\nWidth of detector: %g mm', det_width_mm)
fprintf( '\nWidth of volume: %g mm', vol_width_mm)

%% ASTRA

% Create volume geometry
vol_geom = astra_create_vol_geom( num_voxel_hor, num_voxel_hor, num_voxel_ver);

% Create projection geometry
det_spacing_x = pixel_width_hor_mm / voxel_width_mm;
det_spacing_y = pixel_width_ver_mm / voxel_width_mm;

%det_spacing_x = 1;
%det_spacing_y = det_spacing_x;

%if det_spacing_x > 1
%cfg.option.DetectorSuperSampling  = 1;ceil(det_spacing_x);
%end
%cfg.option.DetectorSuperSampling  = 1;
   
det_row_count = num_pixel_ver;
det_col_count = num_pixel_hor ;
angles = data.fullAngle_rad * ( 0 : num_proj-1 ) / num_proj;
source_origin = source_origin_mm / voxel_width_mm;
origin_det = origin_det_mm / voxel_width_mm;
proj_geom = astra_create_proj_geom('cone', det_spacing_x, det_spacing_y, det_row_count, det_col_count, angles, source_origin, origin_det);
fprintf( '\nRatio: (detector pixel size) / (voxel size): %g, %g', det_spacing_x, det_spacing_y)

% BACKPROJECTION
sino = normat( sino );
[~, rec] = astra_create_backprojection3d_cuda(sino, proj_geom, vol_geom);
rescaleFac = 1 / num_proj / sqrt( num_voxel_hor^2 * num_voxel_ver );
%rec = rescaleFac * lambda * rec;
rec = rescaleFac * rec;

% L-norms
lnorm = @(data) sqrt( sum( data(:).^2 ) );
lnorm_rec = zeros(1,NumIterations);
lnorm_res = zeros(1,NumIterations-1);
lnorm_upd = zeros(1,NumIterations-1);
lnorm_rec(1) = lnorm(rec);

% Print and dispaly
fprintf('\n L2 norms:\n %14s%14s%14s','Solution','Residual','Update')
fprintf('\n %14g', lnorm_rec(1))
subplot(2,2,1,'replace')
zz = round( size( rec, 3) / 2 );
imsc( rec(:,:,zz) )
title('solution')

for nn = 1:(NumIterations-1)

    % FORWARDPROJECTION 
    [~, res] = astra_create_sino3d_cuda(rec, proj_geom, vol_geom);
    res = sino - res;

    % BACKPROJECTION
    [~, recu] = astra_create_backprojection3d_cuda( res, proj_geom, vol_geom);
    recu = rescaleFac * lambda * recu;    
    
    % UPDATE
    rec = rec + recu;
    
    % Dipslay and print
    subplot(2,2,1)
    imsc( rec(:,:,zz) )
    title('solution')
    lnorm_rec(nn+1) = lnorm( rec );
    lnorm_res(nn) = lnorm( res );   
    lnorm_upd(nn) = lnorm( recu );
    fprintf('%14g%14g\n %14g', lnorm_res(nn), lnorm_upd(nn), lnorm_rec(nn+1) )
    % solution
    subplot(2,2,2)    
    plot(lnorm_rec)    
    axis square tight;
    title('l2-norm: solution')
    % residual
    subplot(2,2,3) 
    axis square tight;
    plot(lnorm_res)
    title('l2-norm: residual')
    % update
    subplot(2,2,4)     
    axis square tight;
    plot(lnorm_upd)
    title('l2-norm: update')
end

fprintf('\n')
% Clean up. Note that GPU memory is tied up in the algorithm object,
% and main RAM in the data objects.
astra_clear


function imsc(im)
% imagesc using gray colormap.
imagesc(im)
axis tight equal
colormap(gray)
pause(0.2)
