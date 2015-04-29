function [rec, sino] = astra_fanbeam_test(dataSetNum,showFigs,padding,pixelfilter)

if nargin < 1
    dataSetNum = 4;
end
if nargin < 2
    showFigs = 1;
end
if nargin < 3
    padding = -1;
end
if nargin < 4
    pixelfilter = [0 0 0];
end

num_iterations = 100;
%cfg = astra_struct('SIRT_CUDA');
%cfg = astra_struct('CGLS_CUDA');
cfg = astra_struct('FBP_CUDA');
viewRange = 0.1;
ratio_voxelsToPixels = 1 ;

%% Data

data(4).dir = '/home/jmoosmann/data/gate/';
data(4).filename = '20150317_water_spheres_all_photons_high_act';
data(4).fieldname = 'detector_astra';
data(4).permuteOrder = [2 1 3];
data(4).fullAngle_rad = 2 * pi;
data(4).detectorWidth_mm = 50 ; 
data(4).distance_source_origin_mm = 290; 
data(4).distance_source_det_mm = 300;
data(4).rotation_axis_position = 0;
data(4).rotation_axis_offset = 0; 
data(4).takeNegativeLog = 0;

data(5).dir = [getenv('HOME') '/data/walnut/'];
data(5).filename = 'FullSizeSinograms';
data(5).fieldname = 'sinogram1200';
data(5).permuteOrder = [2 1];
data(5).fullAngle_rad = 2 * pi;
data(5).detectorWidth_mm = 114.8 ; 
data(5).distance_source_origin_mm = 110; 
data(5).distance_source_det_mm = 300;
data(5).rotation_axis_position = 0;
data(5).rotation_axis_offset = -12; 
data(5).takeNegativeLog = 1;


%% Read data
data = data(dataSetNum);
fprintf('\n%s',data.filename)
sino = load( [data.dir data.filename '.mat'] );
sino = single (sino.(data.fieldname) );
sino = permute( sino, data.permuteOrder );

fprintf('\nDimension of input sinograms: ')
fprintf('%u ',size(sino))
sino = sino(:,:,round(size(sino,3)/2));

sino = double( sino );

sino = FilterPixel(sino,pixelfilter(1:2),1,[3 3],pixelfilter(3));

fprintf('\n')
domain(sino)
if data.takeNegativeLog
    if isequal( min( sino(:) ), 0)
        %sino = -log( 1 + sino);
    else
       sino = -log ( sino );
       %sino = -sino;
    end
end
%sino = normat( sino );
domain(sino)
sino = sino -min(sino(:));
sino = single( sino );



% Rotation axis
if data.rotation_axis_offset > 0
    sino = sino(:,1+data.rotation_axis_offset:end);
elseif data.rotation_axis_offset < 0
    sino = sino(:,1:end+data.rotation_axis_offset);
end

num_pixel_0  = size(sino,2);
if isequal( padding, -1 )
    padding = ceil( ( sqrt(2) -1 ) / 2 * num_pixel_0 );
    sino = padarray(sino,[0 padding],'replicate','both');
else
    sino = padarray(sino,[0 padding],'replicate','both');
end

num_proj = size(sino,1);
num_pixel  = size(sino,2);
num_voxel = round( ratio_voxelsToPixels * num_pixel );
% View range for imagesc
viewRange = round( 1 + (num_voxel-1) * [viewRange  1-viewRange] );
viewRange = viewRange(1):viewRange(2);

% Print information
fprintf('\nNumber of detector pixels (1D): %4u',num_pixel_0)
fprintf('\nNumber of detector pixels (1D) after padding: %4u',num_pixel)
fprintf('\nNumber of projection: %4u',num_proj)
fprintf('\nNumber of voxels (1D): %4u',num_voxel)


%% Parameters in physical units
% length in mm 
det_width_mm = data.detectorWidth_mm ; % = num_pixel * 0.05
pixel_width_mm = det_width_mm / num_pixel_0 ;
source_origin_mm = data.distance_source_origin_mm; 
source_det_mm = data.distance_source_det_mm;
origin_det_mm = source_det_mm - source_origin_mm;
vol_width_mm = 2 * source_origin_mm * det_width_mm / ( 2 * source_det_mm + det_width_mm );
voxel_width_mm = vol_width_mm / num_voxel;

%% ASTRA

% Create volume geometry
vol_geom = astra_create_vol_geom( num_voxel);

% Create a data object for the reconstruction
rec_id = astra_mex_data2d( 'create', '-vol', vol_geom);

%% Create projection geometry

% normalize parameters to voxel length for ASTRA projection geometry
%det_width: distance between the centers of two adjacent detector pixels =
%pixel_width. !Do not confuse it with detector width
pixel_width = pixel_width_mm / voxel_width_mm ;
%det_count: number of detector pixels in a single projection
%det_count = num_pixel;
%angles: projection angles in radians
angles = data.fullAngle_rad * ( 0:num_proj-1 ) / num_proj;
%source_origin: distance between the source and the center of rotation
source_origin = source_origin_mm / voxel_width_mm;
%origin_det: distance between the center of rotation and the detector array
origin_det = origin_det_mm / voxel_width_mm;

%proj_geom = astra_create_proj_geom('fanflat', det_width, det_count, angles, source_origin, origin_det);
proj_geom = astra_create_proj_geom( 'fanflat', pixel_width, num_pixel, angles, source_origin, origin_det);

% Create projection object and store projections in it
proj_id = astra_mex_data2d( 'create', '-sino', proj_geom, sino );

% Set up the parameters for a reconstruction algorithm using the GPU
cfg.ReconstructionDataId = rec_id;
cfg.ProjectionDataId = proj_id;

% Create the algorithm object from the configuration structure
alg_id = astra_mex_algorithm('create', cfg);

% Run reconstruction
fprintf('\nAlgorithm: %s',cfg.type)
if isequal(cfg.type(1:3),'FBP')
    tic;
    astra_mex_algorithm('iterate', alg_id, 1);
    fprintf('\nElapsed time for reconstruction: %g\n',toc)
    % Fetch data from ASTRA memory
    rec = astra_mex_data2d('get_single', rec_id);
    if showFigs
        imsc(rec)
    end
else
    fprintf('\n Iteration: \n')    
    for nn = 1:num_iterations
        fprintf(' %u',nn)
        tic;
        astra_mex_algorithm('iterate', alg_id, 1);
        %fprintf(' (%.2g s)',toc)
        % Fetch data from ASTRA memory
        rec = astra_mex_data2d('get_single', rec_id);
        if showFigs
            imsc(rec(viewRange,viewRange));
            axis equal tight off; 
            pause(0.05);
        end
    end
    fprintf('\n')
end

%% Show
if showFigs
    itool(rec)
end

aclear;