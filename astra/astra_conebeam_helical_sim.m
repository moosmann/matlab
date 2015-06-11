function [rec, proj_data] = astra_conebeam_helical_sim( NumProj, AlgType_str, NumIterations)

if nargin < 1
    NumProj = 35;
end
if nargin < 2
    AlgType_str = 'SIRT3D';
end
if nargin < 3
    NumIterations = 50;
end

%% Parameters
% Physical
num_proj = NumProj;
num_pixel_0 = 100;
num_pixel = 2 * round( sqrt(2) / 2 * num_pixel_0);
num_voxel = 1 * num_pixel_0;
%num_voxel = num_pixel;

num_voxel_vec = [num_voxel, num_voxel, num_voxel];
fullAngle_rad = 2 * pi;
% Algorithms
num_iterations = NumIterations;
cfg = astra_struct( [ upper(AlgType_str) '_CUDA' ] );
%cfg.option.MinConstraint = 0;cfg.option.MaxConstraint = 1;
cfg.options.GPUindex = 0;

%% Parameters in physical units
% length in mm 
det_width_mm = 50 ;
pixel_width_mm = det_width_mm / num_pixel_0 ;
source_origin_mm = 50; 
origin_det_mm = 100;
source_det_mm = source_origin_mm + origin_det_mm;
vol_width_mm = 2 * source_origin_mm * det_width_mm / ( 2 * source_det_mm + det_width_mm );
voxel_width_mm = vol_width_mm / num_voxel;
z_range_mm =  vol_width_mm / 3;
z_range =  z_range_mm / voxel_width_mm ;

%% Phantom
data = Ball(num_voxel_vec,[0.35 0.35 0.5],0.1) + Ball(num_voxel_vec,[0.65 0.65 0.5],0.15);
%data = cat( 3, zeros(num_voxel_vec), data, zeros(num_voxel_vec) );

%% Data
% Create volume geometry
vol_geom = astra_create_vol_geom( num_voxel_vec);

% Create projection geometry for circular cone beam
det_spacing_x = pixel_width_mm / voxel_width_mm;
det_spacing_y = det_spacing_x;
det_row_count = num_pixel_0;
det_col_count = num_pixel;
angles = fullAngle_rad * ( 0 : num_proj-1 ) / num_proj;
source_origin = source_origin_mm / voxel_width_mm;
origin_det = origin_det_mm / voxel_width_mm;
proj_geom_circ = astra_create_proj_geom('cone', det_spacing_x, det_spacing_y, det_row_count, det_col_count, angles, source_origin, origin_det);

% Create vectors from circular cone-beam geometry and add an z shift to
% simulate helical cone-beam geometry
vectors = zeros(numel(num_proj), 12);
for nn = 1:num_proj
    
    z = -z_range + 2 * z_range * ( nn - 1 )/ (num_proj-1);
    
    % source
    vectors(nn,1) = sin(proj_geom_circ.ProjectionAngles(nn)) * proj_geom_circ.DistanceOriginSource;
    vectors(nn,2) = -cos(proj_geom_circ.ProjectionAngles(nn)) * proj_geom_circ.DistanceOriginSource;
    vectors(nn,3) = z;
    
    % center of detector
    vectors(nn,4) = -sin(proj_geom_circ.ProjectionAngles(nn)) * proj_geom_circ.DistanceOriginDetector;
    vectors(nn,5) = cos(proj_geom_circ.ProjectionAngles(nn)) * proj_geom_circ.DistanceOriginDetector;
    vectors(nn,6) = z;
    
    % vector from detector pixel (0,0) to (0,1)
    vectors(nn,7) = cos(proj_geom_circ.ProjectionAngles(nn)) * proj_geom_circ.DetectorSpacingX;
    vectors(nn,8) = sin(proj_geom_circ.ProjectionAngles(nn)) * proj_geom_circ.DetectorSpacingX;
    vectors(nn,9) = 0;
    
    % vector from detector pixel (0,0) to (1,0)
    vectors(nn,10) = 0;
    vectors(nn,11) = 0;
    vectors(nn,12) = proj_geom_circ.DetectorSpacingY;
end
proj_geom_heli = astra_create_proj_geom('cone_vec', proj_geom_circ.DetectorRowCount, proj_geom_circ.DetectorColCount, vectors);

% Create projection data from this
tic
[proj_id, proj_data] = astra_create_sino3d_cuda(data, proj_geom_heli, vol_geom);
fprintf('\nTime to create projection data: %g s', toc)

%iplay(proj_data,2,1)

%disp(vectors(1,:)),disp(vectors(1,:) * voxel_width_mm)

%% Reconstruction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a data object for the reconstruction
rec_id = astra_mex_data3d('create', '-vol', vol_geom);

% Set up the parameters for a reconstruction algorithm using the GPU
%cfg = astra_struct('SIRT3D_CUDA');
cfg.ReconstructionDataId = rec_id;
cfg.ProjectionDataId = proj_id;

% Create the algorithm object from the configuration structure
alg_id = astra_mex_algorithm('create', cfg);

% Run algorithm
fprintf('\nAlgorithm: %s \n',cfg.type)
%x = 0.05;x = round( 1 + (num_voxel-1) * [x 1-x] );x = x(1):x(2);
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
        PrintNum(nn,30,' %u')
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
fprintf('\n')