% Resample reconstructions to common voxel size for ML-based segmentation
% using a U-net CNN.
ca
clear all

%% Scans to process
scans = {
    '/asap3/petra3/gpfs/p05/2017/data/11003288/processed/syn136_95L_Mg10Gd_8w/'
    '/asap3/petra3/gpfs/p05/2017/data/11003288/processed/syn137_96L_Mg10Gd_8w/'
    '/asap3/petra3/gpfs/p05/2017/data/11003440/processed/syn34_79R_Mg10Gd_8w/'
    '/asap3/petra3/gpfs/p05/2017/data/11003440/processed/syn44_66L_Mg5Gd_12w/'
    '/asap3/petra3/gpfs/p05/2017/data/11003440/processed/syn96_82L_Mg5Gd_8w/'
    '/asap3/petra3/gpfs/p05/2018/data/11005553/processed/syn034_59R_Mg5Gd_12w/'
    '/asap3/petra3/gpfs/p05/2018/data/11005553/processed/syn035_56L_Mg10Gd_12w/'
    '/asap3/petra3/gpfs/p05/2017/data/11003773/processed/syn107_66R_Mg10Gd_12w/'
    '/asap3/petra3/gpfs/p05/2017/data/11003773/processed/syn101_53L_Mg5Gd_4w/'
    '/asap3/petra3/gpfs/p05/2017/data/11003773/processed/syn104_91L_Mg10Gd_4w/'
    '/asap3/petra3/gpfs/p05/2017/data/11003773/processed/syn105_100AR_Mg10Gd_4w/'
    };

%% Output directory
segpath = '/asap3/petra3/gpfs/p05/2017/data/11003288/processed/segmentation_missing/';
CheckAndMakePath( segpath )

% Loop over beamtimes
for nn = length( scans ):-1:1
    
    proc = scans{nn};
    if proc(end) == '/'
        proc(end) = [];
    end
    [~, name] = fileparts( proc);
    
    beamtime_year = proc(24 + (0:3) );
    beamtime_id = proc(34 + (0:7) );
    
    fprintf( '\n%2u: year %s, beamtime %s, folder %s', nn, beamtime_year, beamtime_id, name);
    
    reco_path = [proc filesep 'reco' ];
    
    %% Find float recos
    d = dir( [ reco_path '/float*' ] );
    disdir =  [ d(:).isdir ];
    floats_namecell = { d( disdir ).name };
    if isempty( floats_namecell )
        continue
    end
    
    %% Choose first float reco
    vol_name = floats_namecell{1};
    vol_path = [reco_path filesep vol_name];
    
    %% Read log file
    log_file = [reco_path '/reco.log'];
    fid = fopen( log_file );
    c = textscan( fid, '%s', 'Delimiter', {'\n', '\r'} );
    c = c{1};
    fclose( fid );
    
    %% Effective pixelsize
    for ll = 1:numel( c )
        t = regexp( c{ll}, 'effective_pixel_size' );
        if t
            cc = textscan( c{ll}, '%*s : %f %*s');
            effective_pixel_size = cc{1};
            break
        end
    end
    
    %% Binning factor
    t = regexp( vol_name, 'bin|Bin' );
    if ~isempty( t )
        bin = str2double( vol_name(t+3) );
    else
        bin = 1;
    end
    
    %% Unique name
    unique_scan_name = sprintf( '%s_%s_%s', beamtime_year, beamtime_id, name );
    
    %% Create parameter struct
    tmp.unique_name = unique_scan_name;
    tmp.bin = bin;
    tmp.effective_pixel_size = effective_pixel_size;
    tmp.full_reco_path = vol_path;
    
    scan_struct(nn) = tmp;
    fprintf( '\n  bin: %u, pixelsize: %f: %s, name: %s', bin, effective_pixel_size, unique_scan_name )
end

%% Save parameter struct %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = [ segpath 'scan_paramter_struct.m'];
save( filename, 'scans' );
fprintf( '\nFinished reading scans.\n' )

%% Stitch and resample %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf( '\n\nResampling and stitching' )

parpath = [ segpath 'cpd/' ];
CheckAndMakePath( parpath )
overviewpath = [parpath 'overview'];
CheckAndMakePath( overviewpath );

for nn = 1:numel( scan_struct )
    
    fprintf( '\n%2u: %s', nn, scan_struct(nn).unique_name )
    unique_name = scan_struct(nn).unique_name;
    full_reco_path = scan_struct(nn).full_reco_path;
    [scan_path, subreco_folder] = fileparts( full_reco_path );
    scan_path = fileparts( scan_path );
    bin = double( scan_struct(nn).bin );
    effective_pixel_size = scan_struct(nn).effective_pixel_size;
    effective_pixel_size_binned = bin * effective_pixel_size;
    
    % Print parameters
    fprintf( '\n scan_path: %s,\n subreco_folder: %s', scan_path, subreco_folder)
    fprintf( '\n  bin: %u', bin )
    fprintf( '\n  effective pixelsizse: %f', effective_pixel_size )
    fprintf( '\n  effective pixelsizse_binned: %f', effective_pixel_size_binned )
    
    % Read volume
    im_struct = dir( [ full_reco_path filesep '*.tif' ] );
    V = imread( [full_reco_path filesep im_struct(1).name] );
    V = zeros( [size(V) numel( im_struct )], 'single' );
    fprintf( '\n Reading.' )
    parfor kk = 1:numel( im_struct )
        filename = [full_reco_path filesep im_struct(kk).name ];
        V(:,:,kk) = imread( filename, 'tif' );
    end
    [x,y,z] = size( V );
    fprintf( ' size: %u %u %u', x, y, z )
    
    %% Resample
    fprintf( ' Resampling.' )
        
    % Original grid
    [x,y,z] = size( V );
    xx = 1:x;
    yy = 1:y;
    zz = 1:z;
    [X,Y,Z] = meshgrid( xx, yy, zz );
    fprintf( '\n Sample point ranges:' )
    fprintf( '\n  x: %f %f', xx(1), xx(end) )
    fprintf( '\n  x: %f %f', yy(1), yy(end) )
    fprintf( '\n  x: %f %f', zz(1), zz(end) )
    
    % Query grid
    k = 5 / effective_pixel_size_binned;
    xq = k:k:x;
    yq = k:k:y;
    zq = k:k:z;
    [Xq,Yq,Zq] = meshgrid( xq, yq, zq );
    fprintf( '\n Query point ranges:' )
    fprintf( '\n  x: %f %f', xq(1), xq(end) )
    fprintf( '\n  x: %f %f', yq(1), yq(end) )
    fprintf( '\n  x: %f %f', zq(1), zq(end) )
    %Vq = interp3( V, Xq, Yq, Zq );
    Vq = interp3( X, Y, Z, V, Xq, Yq, Zq );
    [x,y,z] = size( Vq );
    fprintf( ' resampled size: %u %u %u', x, y, z )
    
    %% Save resampled volume
    fprintf( '\n Saving.' )
    pause( 1 )
    unique_scan_name = scan_struct(nn).unique_name(1:end-2);
    outpath = sprintf( '%s%s', parpath, unique_scan_name );
    CheckAndMakePath( outpath )
    parfor kk = 1:z
        filename = sprintf( '%s/vol_%06u.tif', outpath, kk )
        write32bitTIFfromSingle( filename, Vq(:,:,kk) );
    end
    
    %% Save ortho slices
    xx = round( x / 2 );
    yy = round( y / 2 );
    zz = round( z / 2 );
    % Save ortho slices x
    im = rot90( FilterOutlier( squeeze( Vq(xx,:,:) ), 0.02, '', 0, 0 ), -1 );
    filename = sprintf( '%s/%s_x%06uof%06u.png', overviewpath, unique_scan_name, xx, x );
    imwrite( normat(im), filename );
    % Save ortho slices y
    im = rot90( FilterOutlier( squeeze( Vq(:,yy,:) ), 0.02, '', 0, 0 ),-1 );
    filename = sprintf( '%s/%s_y%06uof%06u.png', overviewpath, unique_scan_name, yy, y );
    imwrite( normat(im), filename );
    % Save ortho slices z
    im = FilterOutlier( squeeze( Vq(:,:,zz) ), 0.02, '', 0, 0 );
    filename = sprintf( '%s/%s_z%06uof%06u.png', overviewpath, unique_scan_name, zz, z );
    imwrite( normat(im), filename );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\nFINISHED.\n' )
