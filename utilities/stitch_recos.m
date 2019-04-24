function s = stitch_recos( scan_path )
% Stich reconstructed volumes using log file information.
%
%
% Written by J. Moosmann

if nargin < 1
    scan_path = '/asap3/petra3/gpfs/p05/2018/data/11004263/processed/syn004_96R_Mg5Gd_8w';
end
if nargin < 2
    reco_subfolder = 'float_rawBin2';
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ca;

% Scans to stitch
scan_struct = dir( [scan_path '*' ] );
num_scans = numel( scan_struct );
% Allocate containing volumes
s(num_scans) = struct;

% loop over scans
for nn = 1:num_scans
    
    % Name and folder
    name = scan_struct(nn).name;
    full_path = [scan_struct(nn).folder filesep name];
    s(nn).full_path = full_path;
    
    %% Reco log
    reco_log = [full_path '/reco/reco.log' ];
    if ~exist( reco_log, 'file' )
        fprintf( '\nReco log not found!\n' )
        break
    end
    fid = fopen( reco_log );
    c = textscan( fid, '%s', 'Delimiter', {'\n', '\r'} );
    c = c{1};
    fclose( fid );
    s(nn).reco_log = reco_log;
    
    %% Projection shape
    for ll = 1:numel( c )
        t = regexp( c{ll}, 'im_shape_raw' );
        if t
            cc = textscan( c{ll}, '%*s : %u %u');
            im_shape_raw = [ cc{1} cc{2} ];
            break
        end
    end
    s(nn).im_shape_raw = im_shape_raw;
    
    %% Binning factor
    for ll = 1:numel( c )
        t = regexp( c{ll}, 'raw_binning_factor' );
        if t
            cc = textscan( c{ll}, '%*s : %u %u');
            bin = [ cc{1} cc{2} ];
            break
        end
    end
    s(nn).bin = bin;
    
    %% ROI
    for ll = 1:numel( c )
        t = regexp( c{ll}, 'raw_roi' );
        if t
            cc = textscan( c{ll}, '%*s : %u %u');
            raw_roi = [ cc{1} cc{2} ];
            break
        end
    end
    s(nn).raw_roi = raw_roi;
    
    %% Effective pixelsize
    for ll = 1:numel( c )
        t = regexp( c{ll}, 'effective_pixel_size' );
        if t
            cc = textscan( c{ll}, '%*s : %f %*s');
            effective_pixel_size = cc{1};
            break
        end
    end
    s(nn).effective_pixel_size = effective_pixel_size;
    
    %% Scan log
    scan_log = [regexprep( full_path, 'processed', 'raw' ) filesep name 'scan.log'];
    if ~exist( scan_log, 'file' )
        fprintf( '\nScan log not found!\n' )
        break
    end
    fid = fopen( scan_log );
    cell_of_lines = textscan( fid, '%s', 'Delimiter', {'\n', '\r'} );
    cell_of_lines = cell_of_lines{1};
    fclose( fid );
    s(nn).scan_log = scan_log;
    
    %% s_stage_z
    for ll = 1:numel( cell_of_lines )
        t = regexp( cell_of_lines{ll}, 's_stage_z' );
        if t
            cc = textscan( cell_of_lines{ll}, '%s', 'Delimiter', {'=',',',' '}, 'CollectOutput', 1, 'MultipleDelimsAsOne', 1);
            s_stage_z = str2double( cc{1}{2} ) * 1000;
            break
        end
    end
    s(nn).s_stage_z = s_stage_z;
    
    %% Print parameters
    fprintf( '\n %s, raw_shape: %u %u, raw_roi: %u %u, bin: %u, z: %f, pixelsizse: %f', name, im_shape_raw, raw_roi, bin, s_stage_z, effective_pixel_size )
    
    %% Read volume
    vol_path = [full_path '/reco/' reco_subfolder];
    im_struct = dir( [ vol_path filesep '*.tif' ] );
    im = imread( [vol_path filesep im_struct(1).name] );
    vol = zeros( [size(im) numel( im_struct )], 'single' );
    fprintf( ' Reading.' )
    parfor kk = 1:numel( im_struct )
        filename = [vol_path filesep im_struct(kk).name ];
        vol(:,:,kk) = imread( filename, 'tif' );
    end
    s(nn).vol = vol;
    
end
%% Absolute positions
fprintf( '\nAbsolute positions: ' )
zoffset = min( s(:).s_stage_z );
for nn = 1:num_scans
    
    % Check raw roi order !!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    bin = s(nn).bin;
    effective_pixel_size_binned = double( bin ) * s(nn).effective_pixel_size;
    
    % First absolute vertical position of image edge
    zval1 = s(nn).s_stage_z - zoffset + double( s(nn).raw_roi(1) - 1) * s(nn).effective_pixel_size;
    s(nn).vert_pos.zval1 = zval1;
    
    % Second absolute vertical position of lower image edge
    zval2 = s(nn).s_stage_z - zoffset + double( s(nn).raw_roi(2) - 1 - bin ) * s(nn).effective_pixel_size;
    s(nn).vert_pos.zval2 = zval2;
    
    % Absolute positions
    v1z = size( s(nn).vol, 3 );
    zz = zval1 + ( 0:v1z - 1 ) * effective_pixel_size_binned;
    s(nn).vert_pos.z = zz;
    
    % Print info
    fprintf( '\n %2u. volume:', nn );
    fprintf( '\n   zval1: %f micron = %f micron = %f pixel', zz(1), zval1, zval1 / effective_pixel_size_binned);
    fprintf( '\n   zval2: %f micron = %f micron = %f pixel', zz(end), zval2, zval2 / effective_pixel_size_binned);
end

%% Plot vertical cuts
yy = round( size( s(1).vol, 2 ) / 2 );
xx = round( size( s(1).vol, 1 ) / 2 );
hxz = figure( 'Name', 'xz');
hyz = figure( 'Name', 'yz');
hroi = figure( 'Name', 'vertical ROI' );
for nn = 1:num_scans
    % XZ
    xz = rot90( squeeze( s(nn).vol(:,yy,:) ) );
    figure(hxz)
    subplot(2,1,nn)
    imsc( xz )
    zz = size( s(nn).vol, 3 );
    ytz = 1:round( zz / 10 ):zz;
    ytzl = s(nn).vert_pos.z(ytz);
    yticks( ytz )
    yticklabels( ytzl )
    axis equal tight
    % YZ
    yz = rot90( squeeze( s(nn).vol(xx,:,:) ) );
    figure(hyz)
    subplot(2,1,nn)
    imsc( yz )
    yticks( ytz )
    yticklabels( ytzl )
    axis equal tight
    
    % ROI regarding noise level
    % Use
    im = xz;
    l = ( std( stdfilt( im(:,1:100) ), 1, 2) ) ;
    l = l(end:-1:1);
    % Filter
    sigma = 15;
    sz = 3*length( l );    % length of gaussFilter vector
    x = linspace(-sz / 2, sz / 2, sz);
    gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
    gaussFilter = gaussFilter / sum (gaussFilter); % normalize
    lf = conv( padarray( l, [length(l) 0], 'symmetric', 'both' ), gaussFilter, 'same');
    lf = lf(length(l) + (1:length(l) ) );
    
    lfx = round( length(l)/4):round( length(l)*3/4);
    t = 1.5 * mean( lf(lfx) );
    % Lower index
    v1 = 1;
    while 1
        if lf(v1) > t
            v1 = v1 + 1;
        else
            break
        end
    end
    s(nn).noisecut1 = v1;
    %Higher index
    v2 = length( l );
    while 1
        if lf(v2) > t
            v2 = v2 - 1;
        else
            break
        end
    end
    s(nn).noisecut2 = v2;
    
    % Plot
    figure( hroi )
    subplot(1,2,nn)
    plot( [l, lf, t*ones(size(l))] )
    title( sprintf( 'xz cut level. left: %u, right: %u', v1, v2 ) )
end
drawnow

%% Overlap
for nn = 1:num_scans - 1
    % Absolute z ranges of subsequent volumes
    zval_vol1 = s(nn).vert_pos.z;
    zval_vol2 = s(nn+1).vert_pos.z;
    % Find overlap within ROI region if possible
    
    if zval_vol1(1) < zval_vol2(1)
        % First volume lower than second volume
        
        % Check if the lower end of the ROI of the second volume is still
        % larger than the upper end of the first volume
        z2roi_pos = s(nn+1).noisecut1;
        z2roi_val = zval_vol2(z2roi_pos);
        if z2roi_val > zval_vol1(end)
            z2roi_pos = zval_vol2(1);
            z2roi_val = zval_vol2(z2roi_pos);
        end
        
        % Check if the upper end of the ROI of the first volume is still
        % larger than the lower end of the second volume
        z1roi_pos = s(nn).noisecut2;
        z1roi_val = zval_vol1(z1roi_pos);
        if z1roi_val < zval_vol2(1)
            z1roi_pos = zval_vol1(end);
            z1roi_val = zval_vol1(z1roi_pos);
        end
    else
        % Second volume lower than first volume
        
        % Check if the upper end of the ROI of the second volume is lower
        % than the lower end of the first volume
        z2roi_pos = s(nn+1).noisecut2;
        z2roi_val = zval_vol2(z2roi_pos);
        if z2roi_val < zval_vol1(1)
            z2roi_pos = zval_vol2(end);
            z2roi_val = zval_vol2(z2roi_pos);
        end
        
        % Check if the lower end of the ROI of the first volume remains
        % below the upper end of the second volume
        z1roi_pos = s(nn).noisecut1;
        z1roi_val = zval_vol1(z1roi_pos);
        if z1roi_val < zval_vol2(1)
            z1roi_pos = zval_vol1(1);
            z1roi_val = zval_vol1(z1roi_pos);
        end
    end
    [~, overlap_roi_pos_vol1] = min( abs( zval_vol1 - z2roi_val ) );
    [~, overlap_roi_pos_vol2] = min( abs( zval_vol2 - z1roi_val ) );
    
    % Stitch level
    zstitch = ( zval_vol1(overlap_roi_pos_vol1) + zval_vol2(overlap_roi_pos_vol2) ) / 2;
    [~, p1stitch] = min( abs( zval_vol1 - zstitch ) );
    [~, p2stitch] = min( abs( zval_vol2 - zstitch ) );
    
    % Print info
    fprintf( '\n %u. overlap:', nn )
    if zval_vol1(1) < zval_vol2(2)
        fprintf( '\n  volume 1:' )
        fprintf( '\n   %-20s %f', 'start', zval_vol1(1) )
        fprintf( '\n   %-20s %f', 'noisecut1', zval_vol1(s(nn).noisecut1) )
        fprintf( '\n   %-20s %f', 'start overlap ROI', zval_vol1(overlap_roi_pos_vol1) )
        fprintf( '\n   %-20s %f', 'noisecut2', zval_vol1(s(nn).noisecut2) )
        fprintf( '\n   %-20s %f', 'end overlap ROI', z1roi_val )
        fprintf( '\n   %-20s %f', 'end', zval_vol1(end) )
        fprintf( '\n  volume 2:' )
        fprintf( '\n   %-20s %f', 'start', zval_vol2(1) )
        fprintf( '\n   %-20s %f', 'start overlap ROI', z2roi_val )
        fprintf( '\n   %-20s %f', 'noisecut1', zval_vol2(s(nn+1).noisecut1) )
        fprintf( '\n   %-20s %f', 'end overlap ROI', zval_vol2(overlap_roi_pos_vol2) )
        fprintf( '\n   %-20s %f', 'noisecut2', zval_vol2(s(nn+1).noisecut2) )
        fprintf( '\n   %-20s %f', 'end', zval_vol2(end) )
    else
        % 
        fprintf( '\n   %u volume. pos: %f, index: %u', nn, zval_vol1(overlap_roi_pos_vol1), overlap_roi_pos_vol1 )
        fprintf( '\n   %u volume. pos: %f, index: %u', nn+1, zval_vol2(overlap_roi_pos_vol2), overlap_roi_pos_vol2 )
    end
    % Print stitch level info
    fprintf( '\n Stitch level:')
    fprintf( '\n   %-30s %f', 'overlap roi start volume 1', zval_vol1(overlap_roi_pos_vol1) )
    fprintf( '\n   %-30s %f', 'central stitch level', zstitch )
    fprintf( '\n   %-30s %f %u', 'stitch level volume 1', zval_vol1(p1stitch), p1stitch )
    fprintf( '\n   %-30s %f %u', 'stitch level volume 2', zval_vol2(p2stitch), p2stitch )
    fprintf( '\n   %-30s %f', 'overlap roi end volume 2', zval_vol2(overlap_roi_pos_vol2) )

    zval_stitch_vol1 = zval_vol1(p1stitch);
    s(nn).stitch.vol1.val = zval_stitch_vol1;
    s(nn).stitch.vol1.pos = p1stitch;
    zval_stitch_vol2 = zval_vol2(p2stitch);
    s(nn).stitch.vol2.val = zval_stitch_vol2;
    s(nn).stitch.vol2.pos = p2stitch;
    
    fprintf( '\n')
    
    %% Plot
    figure('Name', sprintf( 'vertically stitched volume part %u', nn) )
    % XZ
    im1 = rot90( squeeze( s(nn).vol(:,yy,p2stitch:end) ) );
    im2 = rot90( squeeze( s(nn+1).vol(:,yy,1:p1stitch) ) );
    xzstitch = cat( 1, im1, im2 );
    subplot(1,2,1)
    imsc( xzstitch )
    axis equal tight
    title( 'xz plane' )
    drawnow
    
    % YZ
    im1 = rot90( squeeze( s(nn).vol(:,yy,p2stitch:end) ) );
    im2 = rot90( squeeze( s(nn+1).vol(:,yy,1:p1stitch) ) );
    yzstitch = cat( 1, im1, im2 );
    subplot(1,2,2)
    imsc( yzstitch )
    axis equal tight
    title( 'yz plane' )
    drawnow
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\n' )
