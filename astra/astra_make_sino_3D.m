function proj_data = astra_make_sino_3D( vol, angles, detector_spacing )

%% Defaults %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    angles = round( max( [ size( vol, 1), size( vol, 2)] ) );
end
if nargin < 3
    detector_spacing = [1 1];
end

%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isscalar( angles )
    angles = pi * ( 0:angles -1 ) / angles;
end

% Pixel size
DetectorSpacingX = detector_spacing(1);
DetectorSpacingY = detector_spacing(2);

% Volume geometry
det_col_count = size( vol, 1 );
det_row_count = size( vol, 2 );
slices = size( vol, 3 );
vol_geom = astra_create_vol_geom( det_row_count, det_col_count, slices );

% Geometry vector
vectors = zeros( numel(angles), 12);
for nn = 1:numel( angles )
    
    theta = angles( nn );
    
    % source / ray direction
    vectors(nn,1) =  sin( theta );
    vectors(nn,2) = - cos( theta );
    vectors(nn,3) = 0;

    % center of detector
    vectors(nn,4) = - cos( theta );
    vectors(nn,5) = - sin( theta );
    vectors(nn,6) = 0;

    % vector from detector pixel (0,0) to (0,1)
    vectors(nn,7) = cos( theta ) * DetectorSpacingX;
    vectors(nn,8) = sin( theta ) * DetectorSpacingX;
    vectors(nn,9) = 0;

    % vector from detector pixel (0,0) to (1,0)
    vectors(nn,10) = 0;
    vectors(nn,11) = 0;
    vectors(nn,12) = 1 * DetectorSpacingY;

end

% Projection geometry
proj_geom = astra_create_proj_geom( 'parallel3d_vec',  det_row_count, det_col_count, vectors);


[proj_id, proj_data] = astra_create_sino3d_cuda(vol, proj_geom, vol_geom);


% 
% 
% 
% % Volume object and store data
% vol_id = astra_mex_data3d('create', '-vol', vol_geom, vol);
% 
% % Sino object
% sino_id = astra_mex_data3d( 'create', '-proj3d', proj_geom );
% 
% 
% 
% % Configuration struct
% cfg = astra_struct('FP_CUDA');
% cfg.ProjectionDataId = sino_id;
% cfg.VolumeDataId = vol_id;
% cfg.option.GPUindex = 1:gpuDeviceCount;
% 
% % create sinogram
% alg_id = astra_mex_algorithm('create', cfg);
% astra_mex_algorithm('iterate', alg_id);
% % data as single precision, this is more efficient:
% sino = astra_mex_data3d('get_single', v1);
% 
% % Free memory
% astra_mex_data3d('delete', vol_id);
% astra_mex_data3d('delete', sino_id);
% astra_mex_algorithm('delete', alg_id);
 aclear 