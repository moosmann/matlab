function [rec, proj_data, proj_geom_vec] = astra_conebeam_sim(AlgType_str,NumIterations,num_proj,num_pixel_0)

if nargin < 1
    AlgType_str = 'FDK';
end
if nargin < 2
    NumIterations = 20;
end
if nargin < 3
    num_proj = 360;
end
if nargin < 4
    num_pixel_0 = 200;
end



%% Parameters
% Physical
num_pixel = 2 * round( sqrt(2) / 2 * num_pixel_0);
num_voxel = 1 * num_pixel_0;
%num_voxel = num_pixel;


num_voxel_vec = [num_voxel, num_voxel, num_voxel];
fullAngle_rad = 2 * pi;
% Algorithms
num_iterations = NumIterations;
cfg = astra_struct( [ upper(AlgType_str) '_CUDA' ] );
%cfg.option.MinConstraint = 0;
%cfg.option.MaxConstraint = 1;
cfg.options.GPUindex = 0;

%% Parameters in physical units
% length in mm 
det_width_mm = 50 ;
pixel_width_mm = det_width_mm / num_pixel_0 ;
source_origin_mm = 110; 
origin_det_mm = 100;
source_det_mm = source_origin_mm + origin_det_mm;
vol_width_mm = 2 * source_origin_mm * det_width_mm / ( 2 * source_det_mm + det_width_mm );
voxel_width_mm = vol_width_mm / num_voxel;


%% Phantom
data = Ball(num_voxel_vec,[0.35 0.35 0.5],0.1) + Ball(num_voxel_vec,[0.65 0.65 0.5],0.15);

%% Data
% Create volume geometry
vol_geom = astra_create_vol_geom( num_voxel_vec);

% Create projection geometry
det_spacing_x = pixel_width_mm / voxel_width_mm;
det_spacing_y = det_spacing_x;
det_row_count = num_pixel_0;
det_col_count = num_pixel;
angles = fullAngle_rad * ( 0 : num_proj-1 ) / (num_proj);
source_origin = source_origin_mm / voxel_width_mm;
origin_det = origin_det_mm / voxel_width_mm;
proj_geom = astra_create_proj_geom('cone', det_spacing_x, det_spacing_y, det_row_count, det_col_count, angles, source_origin, origin_det);

% Print info
fprintf('\nNumber of detector pixels: %4u',num_pixel_0)
fprintf('\nNumber of projection: %4u',num_proj)
fprintf('\nNumber of voxels: %4u',num_voxel)
fprintf('\nVoxel width: %4g mm',voxel_width_mm)
fprintf('\nPixel width: %4g mm',pixel_width_mm)
fprintf('\nPixel width / voxel width: %g',det_spacing_x)

% 
proj_geom_vec = astra_geom_2vec(proj_geom);

% Create projection data from this
tic
[proj_id, proj_data] = astra_create_sino3d_cuda(data, proj_geom, vol_geom);
fprintf('\nTime to create projection data: %g s', toc)
%fprintf('\n');domain(proj_data(:))

% A*data = proj_data
% < A*data, proj_data > = < proj_data, Aad*proj_data >
[~, rec] = astra_create_backprojection3d_cuda( proj_data, proj_geom, vol_geom);
lnorm = @(data) sqrt( sum( data(:).^2 ) );
projScaleFac = lnorm( data ) / lnorm( rec );
fprintf('\nScale Factor: |x| / |A*Ax|: %g', projScaleFac*num_proj*num_pixel)
projScaleFac2 = lnorm( proj_data ) / sqrt( sum( data(:).*rec(:) )  );
fprintf('\nScale Factor: |y| / <x,A*x>^(1/2): %g', projScaleFac2)

% proj_data_0 = proj_data;
% x = 36;
% proj_data = single( x * normat( proj_data_0 ) );
% proj_data = -log( exp( -proj_data ) ) / x;
%fprintf('\n');domain(proj_data(:))

%astra_mex_data3d( 'set', proj_id, exp(-proj_data)/10000 );
%fprintf('\n');domain(proj_data(:)) 


% Display a single projection image
%figure, imshow(squeeze(proj_data(:,20,:))',[])
%nimplay(proj_data);

% Create a data object for the reconstruction
rec_id = astra_mex_data3d('create', '-vol', vol_geom);

% Set up the parameters for a reconstruction algorithm using the GPU
%cfg = astra_struct('SIRT3D_CUDA');
cfg.ReconstructionDataId = rec_id;

% Create projection object and store projections in it
% Reconstruction is not invariant under a change of the direction of
% rotation
%proj_geom.ProjectionAngles = proj_geom.ProjectionAngles;
%proj_id = astra_mex_data3d('create', '-proj3d', proj_geom );
%astra_mex_data3d( 'set', proj_id, proj_data );

cfg.ProjectionDataId = proj_id;

% Create the algorithm object from the configuration structure
alg_id = astra_mex_algorithm('create', cfg);

% Run algorithm
fprintf('\nAlgorithm: %s \n',cfg.type)
x = 0.05;
x = round( 1 + (num_voxel-1) * [x 1-x] );
x = x(1):x(2);
if isequal(cfg.type(1:3),'FDK')
    tic;
    astra_mex_algorithm('iterate', alg_id, 1);   
    % Fetch data from ASTRA memory
    rec = astra_mex_data3d('get_single', rec_id);
    imsc( rec(:,:,ceil(num_voxel/2)) )
    axis equal tight off;
else
    fprintf('Iteration: \n')    
    tic;
    res = zeros(1,num_iterations);
    for nn = 1:num_iterations
        %fprintf(' %u',nn)
        PrintNum(nn,20,' %u')
        astra_mex_algorithm('iterate', alg_id, 1);
        % Fetch data from ASTRA memory
        rec = astra_mex_data3d('get_single', rec_id);
        %2-norm of the difference between the projection data and the
        %projection of the reconstruction. (The square root of the sum of
        %squares of differences 
        res(nn) = astra_mex_algorithm('get_res_norm', alg_id);
        % Plotting
        subplot(2,2,1)
        imsc( rec( :, :, ceil(num_voxel/2) ) );
        axis equal tight;
        subplot(2,2,2)
        imsc( squeeze( rec( :, ceil(num_voxel/3), : ) ) );
        axis equal tight;
        subplot(2,2,3)
        axis equal tight;
        imsc( squeeze( rec( ceil(num_voxel/3), :, : ) ) );
        axis equal tight;
        subplot(2,2,4)
        plot(res)
        axis square tight;        
        pause(0.05);
    end
    
end
fprintf('\nTime elapsed for reconstruction: %g s\n',toc)
%figure('Name','res: 2-norm of difference between the projection data and the projection of the reconstruction')
%plot(res)


% Clear memory allocated by ASTRA objects
aclear
%itool( rec,  round(num_voxel/2) , cfg.type) 