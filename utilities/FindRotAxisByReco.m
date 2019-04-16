function [rec, m] = FindRotAxisByReco( sino, angles, offset_range )
astra_clear


% Volume geometry
vx = size( sino, 1 );
vy = vx;
vz = size( sino, 3 );

% Detector geometry
num_proj = numel( angles );
det_col_count = size( sino, 1 );
det_row_count = size( sino, 3 );
DetectorSpacingX = 1;
DetectorSpacingY = 1;

% Create ASTRA geometry vector
vectors = zeros( numel( angles ), 12);

rec = zeros( [vx vy numel( offset_range) ] , 'single' );
m = zeros( [2 numel( offset_range) ] , 'single' );

for oo = 1:numel( offset_range )
    
    
    offset = offset_range( oo );
    
    for nn = 1:num_proj
        
        theta = angles(nn);
        
        % ray direction
        vectors(nn,1) = sin( theta );
        vectors(nn,2) = -cos( theta );
        vectors(nn,3) = 0;
        
        % center of detector
        vectors(nn,4) = -offset * cos( theta );
        vectors(nn,5) = -offset * sin( theta );
        vectors(nn,6) = 0;
        
        % vector from detector pixel (0,0) to (0,1)
        vectors(nn,7) = cos( theta ) * DetectorSpacingX;
        vectors(nn,8) = sin( theta ) * DetectorSpacingX;
        vectors(nn,9) = 0;
        
        % vector from detector pixel (0,0) to (1,0)
        vectors(nn,10) = 0;
        vectors(nn,11) = 0;
        vectors(nn,12) = DetectorSpacingY;
        
    end
    
    % ASTRA projection geometry
    proj_geom = astra_create_proj_geom('parallel3d_vec',  det_row_count, det_col_count, vectors);
    
    % ASTRA volume geometry
    vol_geom = astra_create_vol_geom( vy, vx, vz );
    
    % ASTRA volume object
    vol_id = astra_mex_data3d('create', '-vol', vol_geom);
    
    sino_id = astra_create_sino3d_cuda(vol_id, proj_geom, vol_geom);
    
    % Filter sino
    pad_method = 'symmetric';'replicate';0;'none';
    sinof = FilterSinoForBackproj(sino, 1, 'Ram-Lak', pad_method, 'twice');
    astra_mex_data3d('set', sino_id, sinof)
    
    % ASTRA config struct
    cfg = astra_struct('BP3D_CUDA');
    cfg.ProjectionDataId = sino_id;
    cfg.ReconstructionDataId = vol_id;
    astra_mex_algorithm('create', cfg);
    
    % ASTRA create algorithm object from configuration struct
    bp_id = astra_mex_algorithm('create', cfg);
    
    % ASTRA backprojection
    astra_mex_algorithm('iterate', bp_id, 1);
    
    % Fetch data from ASTRA memory
    im = astra_mex_data3d('get_single', vol_id);
    rec(:,:,oo) = im;
    %* pi/(2*length(theta));
    
    mask_rad = 0.95;
    mask_val = 0;
    roi = double( MaskingDisc( im, mask_rad, mask_val) ) * 2^16;

    m(1,oo) = - mean( roi( roi <= 0 ) );
    m(2,oo) = entropy( roi );
    
    
        
end