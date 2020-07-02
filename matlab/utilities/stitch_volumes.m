function [s, vol] = stitch_volumes( scan_path, scan_subfolder, reco_subfolder, crop, save_stitched_volume, stitched_volume_path )
% Stich reconstructed volumes using log file information.
%
% ARGUMENTS
% scan_path : string. absolute path pattern common for all height levels
% scan_subfolder : string. typically 'reco' or 'reco_phase'
% reco_subfolder : string. subfolder to 'reco' to be used for stitching.
%   E.g. 'float_rawBin2'
% crop : bool, crop volumes at stitch level
% save_stitched_volume : bool. save slices of stitchted volume
% stitched_volume_path : string. default: scan_path
%
% RETRURNS
% s : struct containing individual arrays and full information for
%   stitching
% vol : stitched volume array
%
% Written by J. Moosmann

if nargin < 1
    %scan_path = '/asap3/petra3/gpfs/p05/2017/data/11003950/processed/syn22_77L_Mg5Gd_8w';% upward
    scan_path = '/asap3/petra3/gpfs/p05/2018/data/11004263/processed/syn004_96R_Mg5Gd_8w'; % downward
    %scan_path = '/asap3/petra3/gpfs/p05/2018/data/11004263/processed/syn018_35L_PEEK_8w'; %
end
if nargin < 2
    scan_subfolder = 'reco';
end
if nargin < 3
    reco_subfolder = 'float_rawBin2';
end
if nargin < 4
    crop = 1;
end
if nargin < 5
    save_stitched_volume = 1;
end
if nargin < 6
    stitched_volume_path = '';
end

ca;
%% Read parameters and data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\nscan_path: %s', scan_path )
% Scans to stitch
scan_struct = dir( [scan_path '*' ] );
[~, scan_name ] = fileparts( scan_path );
mask = ~strcmp( {scan_struct.name}, scan_name );
scan_struct = scan_struct( mask );
num_scans = numel( scan_struct );
% Allocate containing volumes
s(num_scans) = struct;

% loop over scans
for nn = 1:num_scans
    
    % Name and folder
    name = scan_struct(nn).name;
    full_path = [scan_struct(nn).folder filesep name];
    s(nn).full_path = full_path;
    
    % Reco log
    reco_log = [full_path filesep scan_subfolder filesep 'reco*.log' ];
    rls = dir( reco_log );    
    if ~isempty( rls )
        reco_log = [rls(end).folder filesep rls(end).name];
    else
        fprintf( '\nReco log not found!\n' )
        break
    end
    fid = fopen( reco_log );
    c = textscan( fid, '%s', 'Delimiter', {'\n', '\r'} );
    c = c{1};
    fclose( fid );
    s(nn).reco_log = reco_log;
    
    % Projection shape
    for ll = 1:numel( c )
        %t = regexp( c{ll}, 'im_shape_raw' );
        t = regexp( c{ll}, {'im_shape_raw', 'raw_image_shape', 'raw_image_shape_raw'} );
        %if t
        if sum( [t{:}] )
            cc = textscan( c{ll}, '%*s : %u %u');
            %im_shape_raw = [ cc{1} cc{2} ];
            im_shape_raw = [cc{:}];
            break
        end
    end
    s(nn).im_shape_raw = im_shape_raw;
    
    % Binning factor
    for ll = 1:numel( c )
        % t = regexp( c{ll}, 'raw_binning_factor' );
        t = regexp( c{ll}, {'raw_bin', 'raw_binning_factor', 'raw_binning'} );
        %if t
        if sum( [t{:}] )
            cc = textscan( c{ll}, '%*s : %u %u');
            %bin = [ cc{1} cc{2} ];
            bin = [cc{:}];
            break
        end
    end
    s(nn).bin = bin;
    
    % ROI
    for ll = 1:numel( c )
        t = regexp( c{ll}, 'raw_roi' );
        if t
            cc = textscan( c{ll}, '%*s : %u %u');
            raw_roi = [ cc{1} cc{2} ];
            break
        end
    end
    if ~exist( 'raw_roi', 'var' ) || ( exist( 'raw_roi', 'var' ) && isempty( raw_roi ) )
        raw_roi(1) = 1;
        raw_roi(2) = im_shape_raw(2);
    end
    s(nn).raw_roi = raw_roi;
    
    % Effective pixelsize
    for ll = 1:numel( c )
        t = regexp( c{ll}, 'effective_pixel_size' );
        if t
            cc = textscan( c{ll}, '%*s : %f %*s');
            effective_pixel_size = cc{1};
            break
        end
    end
    s(nn).effective_pixel_size = effective_pixel_size;
    s(nn).effective_pixel_size_binned = effective_pixel_size * double( bin );
    
    % Scan log
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
    
    % s_stage_z
    for ll = 1:numel( cell_of_lines )
        t = regexp( cell_of_lines{ll}, 's_stage_z' );
        if t
            cc = textscan( cell_of_lines{ll}, '%s', 'Delimiter', {'=',',',' '}, 'CollectOutput', 1, 'MultipleDelimsAsOne', 1);
            s_stage_z = str2double( cc{1}{2} ) * 1000;
            break
        end
    end
    s(nn).s_stage_z = s_stage_z;
    
    % Print parameters
    fprintf( '\n%s:', name )
    fprintf( '\n  raw_shape: %u %u', im_shape_raw )
    fprintf( '\n  raw_roi: %u %u', raw_roi )
    fprintf( '\n  bin: %u', bin )
    fprintf( '\n  s_stage_z: %f', s_stage_z )
    fprintf( '\n  effective pixelsizse: %f', effective_pixel_size )
    
    % Read volume
    vol_path = [full_path filesep scan_subfolder filesep reco_subfolder];
    im_struct = dir( [ vol_path filesep '*.tif' ] );
    im = imread( [vol_path filesep im_struct(1).name] );
    vol = zeros( [size(im) numel( im_struct )], 'single' );
    fprintf( '\n  Reading volume:' )
    num_slices = numel( im_struct );
    % Invert or
    filename_cell = {im_struct(end:-1:1).name};
    parfor kk = 1:num_slices
        filename = [vol_path filesep filename_cell{kk} ];
        vol(:,:,kk) = imread( filename, 'tif' );
    end
    s(nn).vol = vol;
    fprintf( '\n size: %u %u %u', size( vol ) )
end

%% Scan order: z-axis pointing downwards !! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if s(1).s_stage_z > s(2).s_stage_z
    upwards = 1;
    dirstr = 'upward';
    stitch_order = num_scans:-1:1;
else
    upwards = 0;
    dirstr = 'downward';
    stitch_order = 1:num_scans;
end
fprintf( '\nScan direction: %s', dirstr );
s(1).upwards = upwards;
s(1).stitch_order = stitch_order;

%% Absolute positions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\nAbsolute positions: ' )
for nn = 1:num_scans
    
    bin = s(nn).bin;
    effective_pixel_size_binned = s(nn).effective_pixel_size_binned;
    
    % First absolute vertical position of image edge
    zval_first = s(nn).s_stage_z + double( s(nn).raw_roi(1) - 1) * s(nn).effective_pixel_size;
    s(nn).zval_first = zval_first;
    
    % Second absolute vertical position of lower image edge
    zval_last = s(nn).s_stage_z + double( s(nn).raw_roi(2) - 1 - bin ) * s(nn).effective_pixel_size;
    s(nn).zval_last = zval_last;
    
    % Absolute positions
    v1z = size( s(nn).vol, 3 );
    zval = zval_first + ( 0:v1z - 1 ) * effective_pixel_size_binned;
    s(nn).zval = zval;
    
    % Print info
    fprintf( '\n %2u. volume. first val: %f (%f)', nn, zval(1), zval_first );
    fprintf( '\n %2u. volume.  last val: %f (%f)', nn, zval(end), zval_last  );
end

%% Remove offset %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if upwards
    zoffset = min( [s(:).zval] );
else
    zoffset = min( [s(:).zval] );
end
fap = figure( 'Name', 'Absolute Positions' );
lgnd = cell( [1 num_scans] );
for nn = 1:num_scans
    s(nn).zval = s(nn).zval - zoffset;
    s(nn).zval_first = s(nn).zval_first - zoffset;
    s(nn).zval_last = s(nn).zval_last - zoffset;
    
    % Plot absolute positions
    figure( fap )
    hold on
    plot( s(nn).zval )
    lgnd{nn} = sprintf( 'volume %u. length: %u', nn, length( zval ) ) ;
end
title( sprintf( '%s scanning (z-axis pointing downwards', dirstr ) )
legend( lgnd )

%% Noise cut level & Plot vertical cuts
yy = round( size( s(1).vol, 2 ) / 2 );
xx = round( size( s(1).vol, 1 ) / 2 );
hxz = figure( 'Name', sprintf( 'xz %s', dirstr ) );
hyz = figure( 'Name', sprintf( 'yz %s', dirstr ) );
hroi = figure( 'Name', 'vertical ROI' );
for nn = 1:num_scans
    % XZ
    xz = ( squeeze( s(nn).vol(:,yy,:) ) );
    figure(hxz)
    subplot(1,num_scans,stitch_order(nn))
    imsc( xz )
    zz = size( s(nn).vol, 3 );
    ztcks = 1:round( zz / 10 ):zz;
    ztcks_label = round( s(nn).zval(ztcks) );
    xtickangle(90)
    xticks( ztcks )
    xticklabels( ztcks_label )
    axis equal tight
    title( sprintf( 'volume %u', nn ) )
    % YZ
    yz = ( squeeze( s(nn).vol(xx,:,:) ) );
    figure(hyz)
    subplot(1,num_scans,stitch_order(nn))
    imsc( yz )
    xticks( ztcks )
    xtickangle(90)
    xticklabels( ztcks_label )
    axis equal tight
    title( sprintf( 'volume %u', nn ) )
    
    % ROI regarding noise level
    im = xz(end-100:end,:);
    l = ( std( stdfilt( im ), 1, 1) ) ;
    l = l(:);
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
    subplot(1,num_scans,nn)
    plot( [l, lf, t*ones(size(l))] )
    title( sprintf( 'xz cut level. left: %u, right: %u', v1, v2 ) )
    axis tight
end
drawnow

%% Overlap
for nn = 1:num_scans - 1
    
    % Absolute z ranges of subsequent volumes
    zval1 = s(nn).zval;
    zval2 = s(nn+1).zval;
    
    % Find overlap (within noisecut ROI if possible)
    if upwards
        
        % First / upper overlap index and value in first volume
        z1ol_pos = s(nn).noisecut1;
        z1ol_val = zval1(z1ol_pos);
        % Check if first/upper overlap edge of first volume is below lower edge of second volume
        if z1ol_val > zval2(end)
            z1ol_pos = zval1(1);
            z1ol_val = zval1(z1ol_pos);
        end
        
        % Second / lower overlap index and value in second volume
        z2ol_pos = s(nn+1).noisecut2;
        z2ol_val = zval2(z2ol_pos);
        % Check if second/lower overlap edge of second volume is above first
        % volume
        if z2ol_val < zval1(1)
            z2ol_pos = zval2(end);
            z2ol_val = zval2(z2ol_pos);
        end
        
    else
        
        % Lower overlap edge in first volume
        z1ol_pos = s(nn).noisecut2;
        z1ol_val = zval1(z1ol_pos);
        % Check if lower overlap edge of first volume is above upper edge
        % of second volume
        if z1ol_val < zval2(1)
            z1ol_pos = zval1(end);
            z1ol_val = zval1(z1ol_pos);
        end
        
        % Upper overlap edge in second volume
        z2ol_pos = s(nn+1).noisecut1;
        z2ol_val = zval2(z2ol_pos);
        % Check if the lower end of the ROI of the second volume is still
        % larger than the upper end of the first volume
        if z2ol_val > zval1(end)
            z2ol_pos = zval2(1);
            z2ol_val = zval2(z2ol_pos);
        end
        
    end
    
    % Stitch level
    zstitch = ( z1ol_val + z2ol_val ) / 2;
    [~, z1stitch_pos] = min( abs( zval1 - zstitch ) );
    [~, z2stitch_pos] = min( abs( zval2 - zstitch ) );
    z1stitch_val = zval1( z1stitch_pos );
    z2stitch_val = zval2( z2stitch_pos );
    
    if upwards
        
        % Check monotonicity
        while z1stitch_val < z2stitch_val
            z1stitch_pos = z1stitch_pos + 1;
            z1stitch_val = zval1( z1stitch_pos );
        end
        
        % Stich indices
        s(nn).zstitch_pos1 = z1stitch_pos;
        if nn == num_scans - 1
            s(nn).zstitch_pos2 = size( s(nn).vol, 3 );
        end
        if nn == 1
            s(nn+1).zstitch_pos1 = 1;
        end
        s(nn+1).zstitch_pos2 = z2stitch_pos;
    else
        
        % Check monotonicity
        while z1stitch_val > z2stitch_val
            z2stitch_pos = z2stitch_pos + 1;
            z2stitch_val = zval2( z2stitch_pos );
        end
        
        % Stich indices
        if nn == 1
            s(nn).zstitch_pos1 = 1;
        end
        s(nn).zstitch_pos2 = z1stitch_pos;
        s(nn+1).zstitch_pos1 = z2stitch_pos;
        if nn == num_scans - 1
            s(nn+1).zstitch_pos2 = size( s(nn+1).vol, 3 );
        end
    end
    
    % Print overlap info
    fprintf( '\n %u. Overlap:', nn )
    if upwards
        fprintf( '\n   %-30s %f', '1. vol lower edge', zval1(end) )
        fprintf( '\n   %-30s %f', '2. vol lower edge',  zval2(end) )
        fprintf( '\n   %-30s %f', '2. vol lower noisecut2', zval2( s(nn+1).noisecut2 ) )
        fprintf( '\n   %-30s %f', 'stitch level', zstitch )
        fprintf( '\n   %-30s %f', '1. vol upper noisecut1', zval1( s(nn).noisecut1 ) )
        fprintf( '\n   %-30s %f', '1. vol upper edge', zval1(1) )
        fprintf( '\n   %-30s %f', '2. vol upper edge', zval2(1) )
    else
        fprintf( '\n   %-30s %f', '1. vol upper edge', zval1(1) )
        fprintf( '\n   %-30s %f', '2. vol upper ddge', zval2(1) )
        fprintf( '\n   %-30s %f', '2. vol noisecut1', zval2( s(nn+1).noisecut1 ) )
        fprintf( '\n   %-30s %f', 'stitch level', zstitch )
        fprintf( '\n   %-30s %f', '1. vol lower noisecut2', zval1( s(nn).noisecut2) )
        fprintf( '\n   %-30s %f', '1. vol lower edge', zval1(end) )
        fprintf( '\n   %-30s %f', '2. vol lower ddge', zval2(end) )
    end
    
    % Print stitch level info
    fprintf( '\n Stitch level:')
    fprintf( '\n   %-30s value: %f', 'stitch level.', zstitch )
    fprintf( '\n   %-30s value: %f, index: %u', 'stitch level volume 1.', z1stitch_val, z1stitch_pos )
    fprintf( '\n   %-30s value: %f, index: %u', 'stitch level volume 2.', z2stitch_val, z2stitch_pos )
    fprintf( '\n')
    
    % Stitch values, positions and ranges
    if upwards
        zz1 = z1stitch_pos:length(zval1);
        zz2= 1:z2stitch_pos;
    else
        zz1 = 1:z1stitch_pos;
        zz2 = z2stitch_pos:length( zval2 );
    end
    
   % Plot
    figure('Name', sprintf( 'vertically stitched volume part %u', nn ) )
    % XZ
    imx1 = squeeze( s(nn).vol(:,yy,zz1) );
    imx2 = squeeze( s(nn+1).vol(:,yy,zz2) );
    % YZ
    imy1 = squeeze( s(nn).vol(xx,:,zz1) );
    imy2 = squeeze( s(nn+1).vol(xx,:,zz2) );
    if upwards
        xzstitch = rot90( cat( 2, imx2, imx1 ), -1 );
        yzstitch = rot90( cat( 2, imy2, imy1 ), -1 );
    else
        xzstitch = rot90( cat( 2, imx1, imx2 ), -1 );
        yzstitch = rot90( cat( 2, imy1, imy2 ), -1 );
    end
    subplot(1,2,1)
    imsc( xzstitch )
    axis equal tight
    title( 'xz plane' )
    subplot(1,2,2)
    imsc( yzstitch )
    axis equal tight
    title( 'yz plane' )
    drawnow
end

%% Crop volumes
if  crop
    fprintf( '\n Crop volumes ' )
    for nn = 1:num_scans
        % volumes
        s(nn).vol = s(nn).vol(:,:,s(nn).zstitch_pos1:s(nn).zstitch_pos2);
        % z values
        s(nn).zval_stitch = s(nn).zval(s(nn).zstitch_pos1:s(nn).zstitch_pos2);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Save stitched volume
if save_stitched_volume
    
    % Output path
    if isempty( stitched_volume_path )
        stitched_volume_path  = scan_path;
        CheckAndMakePath( stitched_volume_path )
    end
    fprintf( '\n Outpath: %s' , stitched_volume_path )
    
    % Loop over volumes
    index_offset = 0;
    for nn = stitch_order
        
        num_proj = size( s(nn).vol, 3 );
        vol = s(nn).vol;
        
        % Loop over projections
        parfor pp = 1:num_proj
            filename = sprintf( '%s/vol_%06u.tif', stitched_volume_path, pp + index_offset);
            im = vol(:,:,pp)
            write32bitTIFfromSingle( filename, im )
        end
        index_offset = index_offset + num_proj;
    end
end

%% Stitch volume
if nargout == 2
    vol = cat(3, s(s(1).stitch_order).vol );
end
