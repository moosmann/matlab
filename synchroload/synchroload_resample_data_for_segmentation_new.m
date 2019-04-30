%ca
clear all

read_scans = 1;

%% Beamtimes
beamtime = { ...
    2016 11001978
    2017 11003950
    2017 11004016
    2017 11003288
    2017 11003440
    2017 11003440
    2018 11004263
    2018 11004936
    2018 11005553};

%% Output directory
segpath = '/asap3/petra3/gpfs/p05/2018/data/11004936/processed/segmentation/';

%% Read scans %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if read_scans
    
    scans = [];
    n_cpd = 0;
    n_load = 0;
    n_other = 0;
    
    % Loop over beamtimes
    for nn = 1:length( beamtime )
        
        % beamtime folder
        beamtime_year = beamtime{nn,1};
        beamtime_id = beamtime{nn,2};
        proc = sprintf( '/asap3/petra3/gpfs/p05/%u/data/%u/processed', beamtime_year, beamtime_id );
        
        % scan folder
        scan_struct = dir( [proc '/*Mg*'] );
        
        fprintf( '\nexp %2u: year %u, beamtime %u, %u folders with Mg', nn, beamtime_year, beamtime_id, length( scan_struct ) );
        
        % Loop over scans
        for mm = 1:length( scan_struct )
            scan_pattern = scan_struct(mm).name;
            scan_path = [ scan_struct(mm).folder filesep scan_pattern ];
            
            if isfolder( scan_path )
                
                reco_path = [scan_path filesep 'reco' ];
                
                if exist( reco_path, 'dir')
                    
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
                    unique_scan_name = sprintf( '%u_%u_%s', beamtime_year, beamtime_id, scan_pattern);
                    
                    %% Create parameter struct
                    tmp.unique_name =unique_scan_name;
                    tmp.bin = bin;
                    tmp.effective_pixel_size = effective_pixel_size;
                    tmp.full_reco_path = vol_path;
                    
                    %% Sort scan type
                    expr = '_\d\d\d$|load|test|push|drill';
                    t = regexp( scan_pattern, expr, 'once');
                    
                    if isempty( t )
                        % CPD
                        n_cpd = n_cpd + 1;
                        scans.cpd(n_cpd) = tmp;
                        tag = 'c';
                    else
                        
                        expr = 'test' ;
                        t = regexp( scan_pattern, expr, 'once' );
                        if isempty( t )
                            % load
                            n_load = n_load + 1;
                            scans.load(n_load) = tmp;
                            tag = 'l';
                        else
                            % rest
                            n_other = n_other + 1;
                            scans.other(n_other) = tmp;
                            tag = 'o';
                        end
                    end
                    
                    fprintf( '\n %2u. bin: %u, pixelsize: %f, tag: %s, name: %s', mm, bin, effective_pixel_size, tag, unique_scan_name )
                    
                end
            end
        end
        
        
        
    end
    
    %% Save parameter struct
    filename = [ segpath 'scan_paramter_struct.m'];
    save( filename, 'scans' );
    fprintf( '\nFinished reading scans.\n' )
end

%% Stitch and resample %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
process.tags = {'l'};
stitch_only = 1;
resample_and_stitch = 1;
if resample_and_stitch
    fprintf( '\n\nResampling and stitching' )
    
    %% Loop over scan types
    for nn = 1:numel( process.tags )
        
        tag = process.tags{nn};
        
        fprintf( '\n Scan type: %s', tag );
        
        switch tag
            case 'c'
                scan_struct_cur = scans.cpd;
                parpath = [ segpath 'cpd/' ];
            case 'l'
                scan_struct_cur = scans.load;
                parpath = [ segpath 'load/' ];
            case 'o'
                scan_struct_cur = scans.other;
                parpath = [ segpath 'other/' ];
        end
        
        CheckAndMakePath( parpath )
        overviewpath = [parpath 'overview'];
        CheckAndMakePath( overviewpath );
        
        %% Loop over scans
        sprev = '';
        found = 0;
        for mm = 1:numel( scan_struct_cur )
            s = scan_struct_cur(mm);
            scan_name = s.unique_name;
            fprintf( '\n %s', s.unique_name )
        end
    end
end

%% Find height scans and stitch
process_height_scans = 0;
process_single_height_scans = 1;
expr = '_a$';
t = regexp( {scan_struct_cur(:).unique_name}, expr, 'once');
for nn = 1:numel( t )
    scan_name = scan_struct_cur(nn).unique_name;
    if ~isempty( t{nn} )
        l = regexp( {scan_struct_cur(:).unique_name}, scan_name(1:end-2) );
        cnt = 0;
        for mm = 1:numel( l )
            if ~isempty( l{mm} )
                cnt = cnt + 1;
            end
        end
        
        if cnt > 1
            if ~process_height_scans
                continue
            end
            fprintf( '\n%s', scan_name )
            fprintf( '\n height levels: %u', cnt )
            full_reco_path = scan_struct_cur(nn).full_reco_path;
            [scan_path, subreco_folder] = fileparts( full_reco_path );
            scan_path = fileparts( scan_path );
            fprintf( '\n scan_path: %s,\n subreco_folder: %s', scan_path, subreco_folder)
            stitch_struct = find_volume_stitch_parameter( scan_path(1:end-2), subreco_folder  );
            effective_pixel_size_binned = stitch_struct(1).effective_pixel_size_binned;
            
            % Resample
            fprintf( ' Resampling.' )
            
            % Concatenate volume
            V = cat(3, stitch_struct(stitch_struct(1).stitch_order).vol );
            zval = cat( 2, stitch_struct(stitch_struct(1).stitch_order).zval_stitch);
            clear stitch_struct
            
            % Original grid
            [x,y,~] = size( V );
            xx = 1:x;
            yy = 1:y;
            zz = zval / effective_pixel_size_binned;
            [X,Y,Z] = meshgrid( xx, yy, zz );
            fprintf( '\n Sample point ranges:' )
            fprintf( '\n  x: %f %f', xx(1), xx(end) )
            fprintf( '\n  x: %f %f', yy(1), yy(end) )
            fprintf( '\n  x: %f %f', zz(1), zz(end) )
            
            % Query grid
            k = 5 / effective_pixel_size_binned;
            xq = k:k:x;
            yq = k:k:y;
            zq = zz(1):k:zz(end);
            [Xq,Yq,Zq] = meshgrid( xq, yq, zq );
            fprintf( '\n Query point ranges:' )
            fprintf( '\n  x: %f %f', xq(1), xq(end) )
            fprintf( '\n  x: %f %f', yq(1), yq(end) )
            fprintf( '\n  x: %f %f', zq(1), zq(end) )
            
            %Vq = interp3( vol, Xq, Yq, Zq );
            Vq = interp3( X, Y, Z, V, Xq, Yq, Zq );
            [x,y,z] =size( Vq );
            
            %% Save resampled volume
            fprintf( ' Saving.' )
            pause( 1 )
            unique_scan_name = scan_struct_cur(nn).unique_name(1:end-2);
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
    end
end
% Single height scans
fprintf( '\n\nSINGLE HEIGHT SCANS' )
expr = '_a$|_b$';
ssc = regexp( {scan_struct_cur(:).unique_name}, expr, 'once');
for nn = 1:numel( ssc )
    unique_name = scan_struct_cur(nn).unique_name;
    if isempty( ssc{nn} )
        
        fprintf( '\n%s', unique_name )
        full_reco_path = scan_struct_cur(nn).full_reco_path;
        bin = double( scan_struct_cur(nn).bin );
        effective_pixel_size = scan_struct_cur(nn).effective_pixel_size;
        effective_pixel_size_binned = bin * effective_pixel_size;
        
        % Print parameters
        fprintf( '\n  bin: %u', bin )
        fprintf( '\n  effective pixelsizse: %f', effective_pixel_size )
        fprintf( '\n  effective pixelsizse_binned: %f', effective_pixel_size_binned )
        
        % Read volume
        im_struct = dir( [ full_reco_path filesep '*.tif' ] );
        vol = imread( [full_reco_path filesep im_struct(1).name] );
        vol = zeros( [size(vol) numel( im_struct )], 'single' );
        fprintf( '\n Reading.' )
        parfor kk = 1:numel( im_struct )
            filename = [full_reco_path filesep im_struct(kk).name ];
            vol(:,:,kk) = imread( filename, 'tif' );
        end
        [x,y,z] = size( vol );
        fprintf( ' size: %u %u %u', x, y, z )
        
        % Resample
        fprintf( '\n Resampling.' )
        k = 5 / effective_pixel_size_binned;
        xq = k:k:x;
        yq = k:k:y;
        zq = k:k:z;
        [Xq,Yq,Zq] = meshgrid( xq, yq, zq );
        Vq = interp3( vol, Xq, Yq, Zq );
        [x,y,z] = size( Vq );
        fprintf( ' resampled size: %u %u %u', x, y, z )
        
        % Save resampled volume
        fprintf( '\n Saving.' )
        pause( 1 )
        outpath = sprintf( '%s%s', parpath, unique_name );
        CheckAndMakePath( outpath )
        fprintf( ' outpath: %s', outpath );
        parfor kk = 1:z
            filename = sprintf( '%s/vol_%06u.tif', outpath, kk )
            write32bitTIFfromSingle( filename, Vq(:,:,kk) );
        end
        
        % Save ortho slices
        xx = round( x / 2 );
        yy = round( y / 2 );
        zz = round( z / 2 );
        % Save ortho slices x
        im = rot90( FilterOutlier( squeeze( Vq(xx,:,:) ), 0.02, '', 0, 0 ), 1 );
        filename = sprintf( '%s/%s_x%06uof%06u.png', overviewpath, unique_name, xx, x );
        imwrite( normat(im), filename );
        % Save ortho slices y
        im = rot90( FilterOutlier( squeeze( Vq(:,yy,:) ), 0.02, '', 0, 0 ),1 );
        filename = sprintf( '%s/%s_y%06uof%06u.png', overviewpath, unique_name, yy, y );
        imwrite( normat(im), filename );
        % Save ortho slices z
        im = FilterOutlier( squeeze( Vq(:,:,zz) ), 0.02, '', 0, 0 );
        filename = sprintf( '%s/%s_z%06uof%06u.png', overviewpath, unique_name, zz, z );
        imwrite( normat(im), filename );
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\nFINISHED.\n' )
