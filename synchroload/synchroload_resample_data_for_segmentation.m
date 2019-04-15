ca
clear all

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
%parpath = '/asap3/petra3/gpfs/p05/2017/data/11003288/processed/segmentation/';
parpath = '/asap3/petra3/gpfs/p05/2018/data/11004936/processed/segmentation/';
overviewpath = [parpath 'overview'];
CheckAndMakePath( parpath )
CheckAndMakePath( overviewpath );

%parpath_load = '/asap3/petra3/gpfs/p05/2018/data/11004936/processed/segmentation_load/';

scans = {};
ss = 0;

% Loop over beamtimes
for nn = length( beamtime ):-1:1
    
    % beamtime folder
    beamtime_year = beamtime{nn,1};
    beamtime_id = beamtime{nn,2};
    proc = sprintf( '/asap3/petra3/gpfs/p05/%u/data/%u/processed', beamtime_year, beamtime_id );
    fprintf( '\n%2u : %s', nn, proc )
    
    % scan folder
    scan_struct = dir( [proc '/*Mg*'] );
    fprintf( ' #folder with Mg: %u', length( scan_struct ) )
    
    % Loop over scans
    for mm = length( scan_struct ):-1:1
        scan_name = scan_struct(mm).name;
        
        % Exclude load tomograms and others using regular expressions
        expr = '_\d\d\d|load|test|push|drill' ;
        t = regexp( scan_name, expr );
        scan_path = [ scan_struct(mm).folder filesep scan_name ];
        if isempty( t ) && isdir( scan_path )
            
            % Store all scans in a cell array
            ss = ss + 1;
            scans{ss} = scan_path;
            
            reco_path = [scan_path filesep 'reco' ];
            
            if exist( reco_path, 'dir')
                
                %% Find float recos
                d = dir( [ reco_path '/float*' ] );
                disdir =  [ d(:).isdir ];
                n = { d( disdir ).name };
                
                %% Choose first float reco
                vol_name = n{1};
                vol_path = [reco_path filesep vol_name];
                
                fprintf( '\n %2u: %s', mm, vol_path )
            
                %% Read log file
                log_file = [reco_path '/reco.log'];
                fid = fopen( log_file );
                c = textscan( fid, '%s', 'Delimiter', {'\n', '\r'} );
                c = c{1};
                fclose( fid );
                
                %% Read effective pixelsize
                for ll = 1:numel( c )
                    t = regexp( c{ll}, 'effective_pixel_size' );
                    if t
                        cc = textscan( c{ll}, '%*s : %f %*s');
                        effective_pixel_size = cc{1};
                        fprintf( ' pixel_size: %f', effective_pixel_size )
                        break
                    end
                end                
                
                %% Check binning factor
                t = regexp( vol_name, 'bin|Bin' );
                if ~isempty( t )
                    bin = str2double( vol_name(t+3) );
                else
                    bin = 1;
                end
                fprintf( ' bin: %u', bin );
                pixel_size = effective_pixel_size * bin;
                
                %% Read volume
                im_struct = dir( [ vol_path filesep '*.tif' ] );
                vol = imread( [vol_path filesep im_struct(1).name] );
                vol = zeros( [size(vol) numel( im_struct )], 'single' );
                fprintf( ' Reading.' )
                parfor kk = 1:numel( im_struct )
                    filename = [vol_path filesep im_struct(kk).name ];
                    vol(:,:,kk) = imread( filename, 'tif' );
                end
                
                %% Resample
                fprintf( ' Resampling.' )
                [x,y,z] =size( vol );
                k = 5 / pixel_size;
                xq = k:k:x;
                yq = k:k:y;
                zq = k:k:z;
                [Xq,Yq,Zq] = meshgrid( xq, yq, zq );
                Vq = interp3( vol, Xq, Yq, Zq );
                [x,y,z] =size( Vq );
                
                %% Save resampled volume
                fprintf( ' Saving.' )
                unique_scan_name = sprintf( '%u_%u_%s', beamtime_year, beamtime_id, scan_name);
                outpath = sprintf( '%s%s', parpath, unique_scan_name );
                CheckAndMakePath( outpath )
                parfor kk = 1:z
                    filename = sprintf( '%s/vol_%04u.tif', outpath, kk )
                    write32bitTIFfromSingle( filename, Vq(:,:,kk) );
                end
                
                %% Save ortho slices                
                xx = round( x / 2 );
                yy = round( y / 2 );
                zz = round( z / 2 );
                % Save ortho slices x
                im = FilterOutlier( squeeze( Vq(xx,:,:) ), 0.02, '', 0, 0 );
                filename = sprintf( '%s/%s_x%04uof%u.png', overviewpath, unique_scan_name, xx, x );
                imwrite( normat(im), filename );
                % Save ortho slices y
                im = FilterOutlier( squeeze( Vq(:,yy,:) ), 0.02, '', 0, 0 );
                filename = sprintf( '%s/%s_y%04uof%u.png', overviewpath, unique_scan_name, yy, y );
                imwrite( normat(im), filename );
                % Save ortho slices z
                im = FilterOutlier( squeeze( Vq(:,:,zz) ), 0.02, '', 0, 0 );
                filename = sprintf( '%s/%s_z%04uof%u.png', overviewpath, unique_scan_name, zz, z );
                imwrite( normat(im), filename );
                
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\nFINISHED.\n' )

if 1

%% Empty type sctructs
scans_cpd = {};
scans_load = {};
scans_other = {};
n_cpd = 0;
n_load = 0;
n_other = 0;

% Loop over beamtimes
for nn = 1:length( beamtime )
    
    % beamtime folder
    beamtime_year = beamtime{nn,1};
    beamtime_id = beamtime{nn,2};
    proc = sprintf( '/asap3/petra3/gpfs/p05/%u/data/%u/processed', beamtime_year, beamtime_id );
    fprintf( '\n%2u : %s', nn, proc )
    
    % scan folder
    scan_struct = dir( [proc '/*Mg*'] );
    fprintf( ' #folder with Mg: %u', length( scan_struct ) )
    
    % Loop over scans
    for mm = length( scan_struct ):-1:1
        scan_name = scan_struct(mm).name;
        scan_path = [ scan_struct(mm).folder filesep scan_name ];
                       
        if isdir( scan_path )
            
            % Store all scans in a cell array
            ss = ss + 1;
            scans{ss} = scan_path;
            
            reco_path = [scan_path filesep 'reco' ];
            
            if exist( reco_path, 'dir')
                
                %% Find float recos
                d = dir( [ reco_path '/float*' ] );
                disdir =  [ d(:).isdir ];
                n = { d( disdir ).name };
                
                %% Choose first float reco
                vol_name = n{1};
                vol_path = [reco_path filesep vol_name];
                
                fprintf( '\n %2u: %s', mm, vol_path )
            
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
                        fprintf( ' pixel_size: %f', effective_pixel_size )
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
                fprintf( ' bin: %u', bin );
                pixel_size = effective_pixel_size * bin;
                
                %% Create parameter struct
                scan_struct.bin = bin;
                scan_struct.effective_pixel_size = effective_pixel_size;
                scan_struct.full_reco_path = vol_path;
                
                %% Sort scan type
                expr = '_\d\d\d|load|test|push|drill' ;
                if isempty( t )
                    n_cpd = n_cpd + 1;
                    scans_cpd(n_cpd) = scan_struct;
                end
                expr = '_\d\d\d|load|push|drill' ;
                t = regexp( scan_name, expr );
                if ~semitpy( t )
                        n_load = n_load + 1;
                        scans_load(n_load) = scan_struct;
                else
                    n_other = n_other + 1;
                    scans_other(n_other) = scan_struct;
                end
            end
        end
    end
end

%% Output directory
parpath_cpd = '/asap3/petra3/gpfs/p05/2018/data/11004936/processed/segmentation_cpd/';
parpath_load = '/asap3/petra3/gpfs/p05/2018/data/11004936/processed/segmentation_load/';
parpath_other = '/asap3/petra3/gpfs/p05/2018/data/11004936/processed/segmentation_load/';

parpath = parpath_cpd;
overviewpath = [parpath 'overview'];
CheckAndMakePath( parpath )
CheckAndMakePath( overviewpath );

% %% Scan folders
% filename = [ parpath 'foldernames.txt'];
% fid = fopen( filename, 'w' );
% for ss = 1:numel( scans )
%     fprintf( fid, '%s\n', scans{ss} );
% end
% fclose( fid );
% filename = [ parpath 'foldernames.m'];
% save( filename, 'scans' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\nFINISHED.\n' )

end