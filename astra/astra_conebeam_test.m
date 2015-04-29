function [rec, sino] = astra_conebeam_test(DataSetNum,AlgType_str,NumIterations,PadHor_PadVer,scaleFactor,VolInit)

if nargin < 1
    DataSetNum = 10;
end
if nargin < 2
    AlgType_str = 'sirt3d';'FDK';
end
if nargin < 3
    NumIterations = 10;
end
if nargin < 4
    PadHor_PadVer = [0 0];
end
if nargin < 5
    scaleFactor = 1;
end
if nargin < 6
    VolInit = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

showFigs = 1;
padding.size_hor = PadHor_PadVer(1);
padding.size_ver = PadHor_PadVer(2);
padding.valueOrMethod = 'replicate';
cfg = astra_struct( [ upper(AlgType_str) '_CUDA' ] );
cfg.option.GPUindex = 0;

num_iterations = NumIterations;
%x = 0.1;

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
fprintf('\nMin/Max of raw projections: [%g, %g]', min( sino(:) ), max( sino(:) ) )

% Normalise
if DataSetNum == 10
    sino = normat( sino );
else
    sino = sino + 1e1;
    sino =  1 - log( (sino) / max( sino(:)) ) ;
    sino = scale_fac * normat( sino );
end
fprintf('\nMin/Max of projections: [%g, %g]', min( sino(:) ), max( sino(:) ) )

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

fprintf( '\nDetector pixels before/after padding [h x v]: %4u x %4u (before),  %4u %4u (after)', num_pixel_hor_0, num_pixel_ver_0, num_pixel_hor, num_pixel_ver )
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

% Create a data object for the reconstruction
rec_id = astra_mex_data3d('create', '-vol', vol_geom);
if VolInit >= 0
    astra_mex_data3d('set', rec_id, VolInit);
end

% Create projection object and store projections in it
proj_id = astra_mex_data3d('create', '-proj3d', proj_geom , 0);
%sino = SubtractMean(sino);
astra_mex_data3d( 'set', proj_id, sino );


%proj_id = astra_mex_data3d('create', '-proj3d', proj_geom, sino );

%astra_mex_data3d('set', proj_id, p);

% Set up the parameters for a reconstruction algorithm using the GPU
cfg.ReconstructionDataId = rec_id;
cfg.ProjectionDataId = proj_id;


if strcmpi( AlgType_str(1:3), 'sir') || strcmpi( AlgType_str(1:3), 'cgl')
% Volume constraints        
    cfg.option.MinConstraint = 0;
    %cfg.option.MaxConstraint = 0.1;
% Super sampling
    if numel( scaleFactor ) >= 2
        cfg.option.DetectorSuperSampling = round( scaleFactor(2) );
        fprintf('\nDetector super sampling: %g', cfg.option.DetectorSuperSampling)
    end     

end

% Create the algorithm object from the configuration structure
alg_id = astra_mex_algorithm('create', cfg);

% Run algorithm
fprintf('\nAlgorithm: %s \n',cfg.type)
if isequal(cfg.type(1:3),'FDK')
    tic;
    astra_mex_algorithm('iterate', alg_id, 1);
    % Fetch data from ASTRA memory
    rec = astra_mex_data3d('get_single', rec_id);
    if showFigs
        subplot(2,2,1)
        imsc( rec(:,:,ceil(num_voxel_ver/2)) )
        subplot(2,2,2)
        imsc( rec(:,:,ceil(num_voxel_ver/2.5)) )
        subplot(2,2,3)
        imsc( squeeze( rec(ceil(num_voxel_hor/2),:,:) ) )
        subplot(2,2,4)
        imsc( squeeze( rec(:,ceil(num_voxel_hor/2),:) ) )
    end
else
    fprintf('Iteration: \n')    
    tic;
    res = zeros(1,num_iterations);
    for nn = 1:num_iterations
        fprintf(' %u',nn)
        astra_mex_algorithm('iterate', alg_id, 1);
        % Fetch data from ASTRA memory
        rec = astra_mex_data3d('get_single', rec_id);
        if showFigs
        %2-norm of the difference between the projection data and the
        %projection of the reconstruction. (The square root of the sum of
        %squares of differences 
        res(nn) = astra_mex_algorithm('get_res_norm', alg_id);
        % Plotting
        subplot(2,2,1)
        imsc( rec( :, :, ceil(num_voxel_ver/2) ) );
        axis equal tight;
        subplot(2,2,2)
        imsc( squeeze( rec( :, ceil(num_voxel_ver/3), : ) ) );
        axis equal tight;
        subplot(2,2,3)
        axis equal tight;
        imsc( squeeze( rec( ceil(num_voxel_hor/3), :, : ) ) );
        axis equal tight;
        subplot(2,2,4)
        plot(res)
        axis square tight;        
        pause(0.05);
        %imsc( rec( x, x, ceil(num_voxel_ver/2) ) );            axis equal tight off;             pause(0.05);
        end
    end   
end
fprintf('\n%s reconstruction. Elapsed time: %g s',cfg.type,toc)

fprintf('\nVolume size: %u x %u x %u', size(rec))
fprintf('\nMin/Max of volume: [%g, %g]', min( rec(:) ), max( rec(:) ) )
fprintf('\nNumber of INFs and NANs: %u, %u', sum( isinf( rec(:) ) ), sum( isnan( rec(:) ) ) ) 
fprintf('\n')

% Clean up. Note that GPU memory is tied up in the algorithm object,
% and main RAM in the data objects.
aclear

function imsc(im)
% imagesc using gray colormap.
imagesc(im)
colormap(gray)
