
%% Parameters
clear all
read_scans = 1;
ignore_case = 0;
process_tags = {'cpd'};
stitch_only = 0;
% Output directory
outpath = '/asap3/petra3/gpfs/p05/2016/data/11001978/processed/segmentation/';
%outpath = '/asap3/petra3/gpfs/p05/2017/data/11003288/processed/segmentation/';

%% Beamtimes
beamtime = { ...
    2016 11001978
    2017 11003950
    2017 11003288
    2017 11003440
    2017 11003773
    2017 11004016
    2018 11004263
    2018 11004936
    2018 11005553
    2019 11005842
    2019 11006704
    2019 11006991
    };

samples = synchroload_samples;
all_counter = 0;

%% Loop over cpd and load scans
for exp_type = { 'cpd', 'load'}
    et = exp_type{1};
    et_counter = 0;
    
    %% Loop over implant types
    for imp_type = { 'Mg5Gd', 'Mg10Gd', 'Ti', 'PEEK' }
        it = imp_type{1};
        it_counter = 0;
        fprintf( '\n\n%s %s', et, it )
        
        name_pattern = samples.(et).(it);
        % add underscore for unique match using regexp
        for np = 1:numel( name_pattern )
            name_pattern{np} = [ '_' name_pattern{np}];
            sample_name = [ 'c' name_pattern{np} ];
            samples.control.(et).(it).( sample_name ) = 0;
        end
        
        %% Loop over beamtimes
        for bt = 1:size( beamtime, 1 )
            
            % beamtime folder
            if size( beamtime, 1 ) == 1
                beamtime_year = beamtime{1};
                beamtime_id = beamtime{2};
            else
                beamtime_year = beamtime{bt,1};
                beamtime_id = beamtime{bt,2};
            end
            proc = sprintf( '/asap3/petra3/gpfs/p05/%u/data/%u/processed', beamtime_year, beamtime_id );
            
            %%
            % scan folder
            dir_struct = dir( proc );
            match2 = contains( {dir_struct.name}, name_pattern, 'IgnoreCase', ignore_case );
            match = zeros( [1, numel(dir_struct)], 'logical' );
            sample_arr = cell( [1, numel(dir_struct)] );
            for mm = 1:numel( match )
                for np = 1:numel(name_pattern )
                    name_pattern_match = contains( dir_struct(mm).name, name_pattern{np}, 'IgnoreCase', ignore_case );
                    sample_name = [ 'c' name_pattern{np} ];
                    if name_pattern_match
                        match(mm) =  1;
                        samples.control.(et).(it).( sample_name ) = 1;
                        sample_arr{mm} = name_pattern{np};
                    end
                end
            end
            if ~isempty( match ) && ~isempty( match2 ) && ~isequal( match, match2 )
                error( '\n Something is wrong' )
            end
            %%
            
            scan_struct = dir_struct(match);
            fprintf( '\nexp %2u %s %s %u %u %u: %s found', bt, et, it, beamtime_year, beamtime_id, length( scan_struct ) );
            sample_arr = sample_arr(match);
            
            % Loop over scans
            for mm = 1:length( scan_struct )
                bin = [];
                effective_pixel_size = [];
                camera = [];
                im_shape_raw = [];
                folder_name = scan_struct(mm).name;
                scan_path = [ scan_struct(mm).folder filesep folder_name ];
                sample = sample_arr{mm};
                
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
                        subfolder_reco = floats_namecell{1};
                        data_path = [reco_path filesep subfolder_reco];
                        
                        %% Read log file
                        dir_struct = dir( [reco_path '/reco*.log'] );
                        log_file = [ reco_path  filesep  dir_struct(1).name ];
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
                        %% Energy
                        for ll = 1:numel( c )
                            t = regexp( c{ll}, 'energy' );
                            if t
                                cc = textscan( c{ll}, '%*s : %f %*s');
                                energy = cc{1};
                                break
                            end
                        end
                        %% Camera
                        for ll = 1:numel( c )
                            t = regexp( c{ll}, 'camera' );
                            if t
                                cc = textscan( c{ll}, '%*s : %s');
                                camera = cc{1}{1};
                                break
                            end
                        end
                        %% im shape raw
                        for ll = 1:numel( c )
                            t = regexp( c{ll}, {'im_shape_raw', 'raw_image_shape', 'raw_image_shape_raw'} );
                            %disp( c{ll} )
                            if sum( [t{:}] )
                                cc = textscan( c{ll}, '%*s : %u %u');
                                im_shape_raw = [cc{:}];
                                break
                            end
                        end
                        %% Binning factor
                        for ll = 1:numel( c )
                            t = regexp( c{ll}, {'raw_bin', 'raw_binning_factor', 'raw_binning'} );
                            if sum( [t{:}] )
                                cc = textscan( c{ll}, '%*s : %u %u');
                                bin = [cc{:}];
                                break
                            end
                        end
                        if ~isscalar( bin )                            
                            t = regexp( subfolder_reco, 'bin|Bin' );
                            if ~isempty( t )
                                bin = str2double( subfolder_reco(t+3) );
                            else
                                bin = 1;
                            end
                        end
                        %% Unique and descriptive name
                        name = sprintf( '%s_%s_E%05.0feV_cam%ux%u_pix%04.0fnm%s_%u_%u_%s', et, it, energy, im_shape_raw, 1000*effective_pixel_size, sample,  beamtime_year, beamtime_id, folder_name);
                        %% Create parameter struct
                        tmp.type = et;
                        tmp.material = it;
                        tmp.energy = energy;
                        tmp.im_shape_raw = im_shape_raw;
                        tmp.effective_pixel_size = effective_pixel_size;
                        tmp.sample = sample(2:end);
                        tmp.beamtime_year = beamtime_year;
                        tmp.beamtime_id = beamtime_id;
                        tmp.folder_name = folder_name;
                        tmp.camera = camera;
                        tmp.bin = bin;
                        tmp.name = name;
                        tmp.data_path = data_path;

                        % Assign to groups
                        it_counter = it_counter + 1;
                        scans.(et).(it)(it_counter) = tmp;
                        et_counter = et_counter + 1;
                        scans.(et).all(et_counter) = tmp;
                        all_counter = all_counter + 1;
                        scans.all(all_counter) = tmp;
                        fprintf( '\n %2u. bin: %u, effpix: %f, cam: %s, shape: %u %u, name: %s', mm, bin, effective_pixel_size, camera, im_shape_raw, name )
                    end
                end
            end
        end
    end % imp_type
end % exp_type
%% Save parameter struct
CheckAndMakePath( outpath )
mfile_path = [outpath 'matlab/'];
CheckAndMakePath( mfile_path )
filename = [ mfile_path 'scans.m'];
save( filename, 'scans' );
filename = [ mfile_path 'samples.m'];
save( filename, 'samples' );
fprintf( '\nFinished reading scans.\n' )


%         % Write CVS file
%         csvfilename=[out_path filesep struct_scans(nn).name '_waiting.csv'];
%         cHeader = {'wait_time' 'wait_force' 'wait_pusher'}; %dummy header
%         commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
%         commaHeader = commaHeader(:)';
%         textHeader = cell2mat(commaHeader); %cHeader in text with commas
%         
%         %write header to file
%         fid = fopen(csvfilename,'w');
%         fprintf(fid,'%s\n',textHeader);
%         fclose(fid);
%         
%         datvalues(:,1)=(wait_time-wait_time(1)) / 60 / 60;
%         datvalues(:,2)=wait_force;
%         datvalues(:,3)=wait_pusher;
%         
%         %write data to end of file
%         dlmwrite(csvfilename,datvalues,'-append');
        

%% Data sets to process %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\n\nData sets to process' )
% Loop over scan types
for nn = 1:numel( process_tags )
    tag = process_tags{nn};
    fprintf( '\n Scan type: %s', tag );
    scans_proc = scans.(tag).all;
    % Loop over scans
    for mm = 1:numel( scans_proc )
        s = scans_proc(mm);
        scan_name = s.name;
        fprintf( '\n %s', s.name )
    end
end
overviewpath = [outpath 'overview'];
CheckAndMakePath( overviewpath );
parpath = [outpath 'scans/'];
CheckAndMakePath( parpath );
        
%% Find height scans and stitch
if 1
    sr_counter = 0;
    expr = '_b$';
    t = regexp( {scans_proc(:).name}, expr, 'once');
    for nn = 127:numel( t )
        scan = scans_proc(nn);
        scan_name = scan.name;
        if ~isempty( t{nn} )
            l = regexp( {scans_proc(:).name}, scan_name(1:end-2) );
            cnt = 0;
            for mm = 1:numel( l )
                if ~isempty( l{mm} )
                    cnt = cnt + 1;
                end
            end
            if cnt > 1
                fprintf( '\n%s', scan_name )
                fprintf( '\n height levels: %u', cnt )
                data_path = scan.data_path;
                [scan_path, subreco_folder] = fileparts( data_path );
                scan_path = fileparts( scan_path );
                fprintf( '\n scan_path: %s,\n subreco_folder: %s', scan_path, subreco_folder)
                stitch_struct = stitch_volumes( scan_path(1:end-2), 'reco', subreco_folder  );
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
                clear X Y Z V Xq Yq Zq;
                [x,y,z] = size( Vq );
                
                % Check volume
                if sum( isinf( Vq(:) ) ) || sum( isnan( Vq(:) ) )
                    error( 'Check volume' )
                end
                vmin = min( Vq(:) );
                vmax = max( Vq(:) );
                vmean = mean( Vq(:) );
                vstd = std( Vq(:) );
                scan.vmin = vmin;
                scan.vmax = vmax;
                scan.vmean = vmean;
                scan.vstd = vstd;
                fprintf( '\n [min, max, max-min, mean, std] = [%f, %f, %f, %f, %f] ', vmin, vmax, vmax - vmin, vmean, vstd )
                [N,edges] = histcounts( Vq(:), 1024 );
                scan.histcounts.N = N;
                scan.histcounts.edges = edges;
                unique_scan_name = scan.name(1:end-2);
                sr_counter = sr_counter + 1;
                scans_resampled( sr_counter ) = scan;
                scans_resampled( sr_counter ).name = unique_scan_name;
                
                %% Save resampled volume
                fprintf( ' Saving.' )
                pause( 1 )
                save_path = sprintf( '%s%s', parpath, unique_scan_name );
                CheckAndMakePath( save_path )
                parfor kk = 1:z
                    filename = sprintf( '%s/vol_%06u.tif', save_path, kk )
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
                clear Vq
            end
        end
    end
    % Single height scans
    fprintf( '\n\nSINGLE HEIGHT SCANS' )
    expr = '_a$|_b|_c$';
    ssc = regexp( {scans_proc(:).name}, expr, 'once');
    for nn = 1:numel( ssc )
        name = scans_proc(nn).name;
        if isempty( ssc{nn} )
            
            fprintf( '\n%s', name )
            data_path = scans_proc(nn).data_path;
            bin = double( scans_proc(nn).bin );
            effective_pixel_size = scans_proc(nn).effective_pixel_size;
            effective_pixel_size_binned = bin * effective_pixel_size;
            
            % Print parameters
            fprintf( '\n  bin: %u', bin )
            fprintf( '\n  effective pixelsizse: %f', effective_pixel_size )
            fprintf( '\n  effective pixelsizse_binned: %f', effective_pixel_size_binned )
            
            % Read volume
            im_struct = dir( [ data_path filesep '*.tif' ] );
            vol = imread( [data_path filesep im_struct(1).name] );
            vol = zeros( [size(vol) numel( im_struct )], 'single' );
            fprintf( '\n Reading.' )
            parfor kk = 1:numel( im_struct )
                filename = [data_path filesep im_struct(kk).name ];
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

            % Check volume
            if sum( isinf( Vq(:) ) ) || sum( isnan( Vq(:) ) )
                error( 'Check volume' )
            end
            vmin = min( Vq(:) );
            vmax = max( Vq(:) );
            vmean = mean( Vq(:) );
            vstd = std( Vq(:) );
            scan.vmin = vmin;
            scan.vmax = vmax;
            scan.vmean = vmean;
            scan.vstd = vstd;
            fprintf( '\n [min, max, max-min, mean, std] = [%f, %f, %f, %f, %f] ', vmin, vmax, vmax - vmin, vmean, vstd )
            [N,edges] = histcounts( Vq(:), 1024 );
            scan.histcounts.N = N;
            scan.histcounts.edges = edges;
            unique_scan_name = scan(nn).name(1:end-2);
            sr_counter = sr_counter + 1;
            scans_resampled( sr_counter ) = scan;
            scans_resampled( sr_counter ).name = unique_scan_name;
            
            %% Save resampled volume
            fprintf( '\n Saving.' )
            pause( 1 )
            save_path = sprintf( '%s%s', parpath, unique_scan_name );
            CheckAndMakePath( save_path )
            fprintf( '\n save_path: %s', save_path );
            parfor kk = 1:z
                filename = sprintf( '%s/vol_%06u.tif', save_path, kk );
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
filename = [ mfile_path 'scans_resampled.m'];
save( filename, 'scans_resampled' );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\nFINISHED.\n' )
