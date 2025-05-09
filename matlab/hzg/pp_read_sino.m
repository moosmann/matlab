%% Read sinogram script %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if par.read_sino
    fprintf('\nRead sino:')
    t = toc;
    % Sino path
    if ~isempty( par.read_sino_folder )
        write.sino_path = [scan_path par.read_sino_folder filesep];
    end
    sino_path = write.sino_path;
    
    % Sinogram file names
    data_struct = dir( [sino_path filesep '*.tif'] );
    if isempty( data_struct )
        error('\n Sinograms not found.')
    end
    sino_names = {data_struct.name};    
    sino_names_mat = NameCellToMat( sino_names );
    num_sino_found = numel(sino_names);
    
    % Parameters
    filename = sprintf('%s%s', sino_path, sino_names_mat( round(num_sino_found / 2 ),:));
    par.tifftrafo = 0;
    [sino, par.tiff_info] = read_image( filename, par );
    read_sino_trafo = par.read_sino_trafo;
    sino = read_sino_trafo( sino );
    
    switch numel(par.sino_roi)
        case 0
            x0 = 1;
            x1 = size( sino, 2 );
            par.im_shape_cropbin1 = size( sino, 2 );
        case 1
            x0 = par.sino_roi(1);
            x1 = size(sino,2) - par.sino_roi-1;
            par.im_shape_cropbin1 = x1 - x0 + 1;
        case 2
            x0 = par.sino_roi(1);
            x1 = par.sino_roi(2);
            par.im_shape_cropbin1 = x1 - x0 + 1;
    end
    num_proj_read = size( sino, 1);
    num_proj_used = size( sino, 1 );
    
    read_sino_range = par.read_sino_range;
    if isempty(read_sino_range)
        read_sino_range = 1;
    end
    switch numel(read_sino_range)
        case 1
            if read_sino_range >= 0 && read_sino_range < 1
                read_sino_range = floor( (num_sino_found - 1) * read_sino_range) + 1;
                fprintf('\n sinogram slice: %u of %u', read_sino_range, num_sino_found)
            elseif read_sino_range >= 1 && read_sino_range < 10
                read_sino_range = 1:read_sino_range:num_sino_found;
            end
        case 2
            if read_sino_range(1) >= 0 && read_sino_range(1) < 1
                read_sino_range(1) = floor( (num_sino_found - 1) * read_sino_range(1) + 1 );
            end
            if read_sino_range(2) >= 0 && read_sino_range(2) < 1
                read_sino_range(2) = floor( (num_sino_found - 1) * read_sino_range(2) + 1 );
            end
            if read_sino_range(2) < 0 && mod(read_sino_range(2),1) == 0
                read_sino_range(2) = num_sino_found + read_sino_range(2);
            end
            fprintf('\n sinogram rang: [%u %u] of %u', read_sino_range(1), read_sino_range(2), num_sino_found)
            read_sino_range = read_sino_range(1):read_sino_range(2);
    end
    sino_names_mat = sino_names_mat(read_sino_range,:);
    num_sino_used = size(sino_names_mat,1);
    im_shape_binned2 = num_sino_used;
    par.num_sino_used = num_sino_used;
    par.num_sino_found = num_sino_found;
    par.im_shape_binned2 = num_sino_used;
    
    %% Read parameters from reconlog.txt
    fn = [scan_path filesep 'reconlog.txt'];
    if exist(fn, 'File')
        fprintf('\nReading reconlog.txt:')
        fid = fopen(fn);
        c = textscan(fid, '%s : %s');
        fclose(fid);
        
        % raw_bin
        b = contains(c{1},'rawbin');
        if sum(b)
            rawbin = c{2}(b);
             rawbin = str2double(rawbin{1});
            fprintf('\n raw bin: %.1f',rawbin);
        end
        par.raw_bin = rawbin;
        
        % angle increment
        b = contains(c{1},'scan_angularstep');
        if sum(b)
            angularstep = c{2}(b);
            angularstep = pi / 180 *str2double(angularstep{1});
            fprintf('\n angle increment: %f rad = %f deg', angularstep, angularstep * 180 / pi)
        end
        
        % number of projections
        b = contains(c{1},'scan_n_angle');
        if sum(b)
            n_angle = c{2}(b);
            n_angle = str2double(n_angle{1});
            fprintf('\n number of angles: %u', n_angle)
            angles = (0:n_angle - 1) * angularstep;
            ar = n_angle * angularstep;
            if isempty(tomo.rot_angle_full_range)
                tomo.rot_angle_full_range = ar;
            end
        end
        
        %scan  pixel size
        b = contains(c{1},'scan_pixelsize');
        if sum(b)
            pixelsize = c{2}(b);
            pixelsize = str2double(pixelsize{1});            
            fprintf('\n scan pixelsize: %f micron', pixelsize * 1000)
            c{1}(b) = [];
            c{2}(b) = [];
        end
        
        % pixel size
        b = contains(c{1},'pixelsize');
        if sum(b)
            pixelsize = str2double(c{2}{b});
            fprintf('\n pixelsize: %f micron', pixelsize * 1000)
            eff_pixel_size = rawbin * pixelsize / 1000;
        end
        
        % reco mode
        b = contains(c{1},'reco_mode');
        if sum(b)
            reco_mode = str2double(c{2}{b});
            fprintf('\n reco mode: %f', reco_mode)
            ar =  reco_mode * pi / 180;
        end
        
        if isempty(tomo.rot_angle_full_range)
            tomo.rot_angle_full_range = ar;
        end
        
        if isempty(eff_pixel_size)
            eff_pixel_size = raw_bin * pixelsize;
        end
        par.eff_pixel_size_binned = eff_pixel_size;
        par.eff_pixel_size = eff_pixel_size / rawbin;
        phase_retrieval.eff_pixel_size_binned = par.eff_pixel_size_binned;
        phase_retrieval.eff_pixel_size = par.eff_pixel_size;
        
    end
    ar = tomo.rot_angle_full_range;
    fprintf('\n angular range: %f rad = %f * pi rad = %f deg', ar, ar / pi, ar * 180 / pi)
    
    % Angles
    if  ~exist( 'angles', 'var' )
        if isempty( tomo.rot_angle_full_range )
            cprintf( 'Red', '\nEnter full angle of rotation (including one additional increment) or vector of angles, in radians: ' );
            tomo.rot_angle_full_range = input( '' );
        end
        if isscalar( tomo.rot_angle_full_range )
            angles = tomo.rot_angle_full_range * (0:num_proj_used - 1) / num_proj_used;
        else
            angles = tomo.rot_angle_full_range;
        end
        if length( angles ) ~= num_proj_used
            error( 'Number of angles (%u) entered not consistent with sinogram (%u) read.', numel( angles), num_proj_read )
        end
    end
    if tomo.rot_angle_full_range > 90
        error('rotation angle full range, %f, is most likely given in degrees instead of rad', tomo.rot_angle_full_range)
    end
    
    if isempty(par.sample_detector_distance)
        d = dir([par.nexus_path filesep '*.h5']);
        nexuslog_name = [d.folder filesep d.name];
        par.sample_detector_distance = h5read( nexuslog_name, '/entry/scan/setup/o_ccd_dist');
        par.sample_detector_distance = par.sample_detector_distance / 1000;
        fprintf('\n sample dector distance: %.3f m', par.sample_detector_distance )
    end
    
    if isempty(par.energy)
        d = dir([par.nexus_path filesep '*.h5']);
        nexuslog_name = [d.folder filesep d.name];
        %% TODO: Fix for P05 and different entry names
        par.energy= h5read( nexuslog_name, '/entry/scan/setup/p07_energy');
        fprintf('\n energy: %.1f keV', par.energy / 1e3 )
        
    end
    par.wave_length = E_to_lambda( par.energy);
    fprintf('\n wave length lambda: %.2g angstrom', par.wave_length * 10^10 )
    par.fresnel_number = par.eff_pixel_size_binned^2 / (par.wave_length * par.sample_detector_distance);
    fprintf('\n fresnel number =  pixelsize^2 / lambda / z: %f',  par.fresnel_number )
    
    %% Preallocation
    fprintf( '\nRead sinograms.')    
    proj = zeros( par.im_shape_cropbin1, num_sino_used, num_proj_read, 'single');
    fprintf( '\n sinograms found: %u', num_sino_found);
    fprintf( '\n sinograms used: %u', num_sino_used);
    fprintf( '\n allocated memory: %.2f GiB.', Bytes( proj, 3 ) )
    fprintf( '\n shape: [num_pixel num_slice num_proj] = [%u %u %u]', par.im_shape_cropbin1, num_sino_used, num_proj_read )
    
    %% Read sinogram
    im_shape_cropbin1 = par.im_shape_cropbin1;
    im_shape_binned1 = im_shape_cropbin1;
    filter_sino = par.filter_sino;    
    gpus = par.gpu_index;
    num_gpu = numel(gpus);        
    if par.use_gpu_in_parfor
        % Limit by GPU pixel filter: 2 * single + 2 * uint16 + 1 * logical
        gpu_mem_requ_per_im = prod( size(sino) + 2 * pixel_filter_radius ) * (4 + 4 + 2 + 2 + 1 );
        poolsize_max_gpu = floor( par.poolsize_gpu_limit_factor * min( mem_avail_gpu ) / gpu_mem_requ_per_im );
        %% TODO: Check max poolsize
        %poolsize_max_gpu = max([poolsize_max_gpu,numel(par.gpu_index)]);
        poolsize_max_gpu = poolsize_max_gpu * numel(par.gpu_index);
        fprintf(' \n estimated GPU memory required per image for pixel filtering : %g MiB', gpu_mem_requ_per_im / 1024^2 )
        fprintf(' \n GPU poolsize limit factor : %g', par.poolsize_gpu_limit_factor )
        fprintf(' \n GPU memory induced maximum poolsize : %u', poolsize_max_gpu )
    else
        poolsize_max_gpu = par.poolsize;
    end
    strong_abs_thresh = par.strong_abs_thresh;
    fprintf('\n strong absorption threshold: %g',strong_abs_thresh)
    if ndims(sino)==1
        filename = sprintf('%s%s', sino_path, sino_names_mat);
        sino = read_image(filename);
        sino = read_sino_trafo(sino);
        sino = sino(x0:x1,:);
        if filter_sino
            sino = FilterPixel(sino,pixel_filter_sino);
        end
        if strong_abs_thresh < 1
            m = sino < strong_abs_thresh;
            sino(m) = mean2(sino(m));
        end
        proj(:, 1,:) = reshape( sino, [im_shape_cropbin1, 1, num_proj_read] );
    else
        parfor (nn = 1:num_sino_used,poolsize_max_gpu)
            filename = sprintf('%s%s', sino_path, sino_names_mat(nn,:));
            sino = read_image(filename);
            sino = read_sino_trafo(sino);
            sino = sino(x0:x1,:);
            if filter_sino
                gpu_index = gpus(mod(nn, num_gpu) + 1);
                sino = FilterPixel(sino,pixel_filter_sino,gpu_index);
            end
            if strong_abs_thresh < 1
                m = sino < strong_abs_thresh;
                sino(m) = mean2(sino(m));
            end
            proj(:, nn,:) = reshape( sino, [im_shape_cropbin1, 1, num_proj_read] );
        end
    end
    fprintf( '\n duration %.1f s (%.2f min)', toc - t, ( toc - t ) / 60 )

    %% Plot data
    if par.visual_output
        t = toc;
        f = figure( 'Name', 'projection and sinogram', 'WindowState', par.window_state );
        fprintf('\nDisplay and save figure: %s',f.Name)
        subplot(1,2,1)
        imsc1( proj(:,:,1))
        title(sprintf('intensity: proj(:,:,1)'))
        colorbar
        axis equal tight
        
        subplot(1,2,2)
        nn = round( size( proj, 2) / 2);
        sino = squeeze( proj(:,nn,:) );
        %sino = FilterOutlier( sino );
        imsc1( sino )
        title(sprintf('sinogram: proj(:,%u,:)', nn))
        colorbar
        axis equal tight
        
        drawnow
        
        CheckAndMakePath( write.fig_path )
        saveas( f, sprintf( '%s%s.png', write.fig_path, regexprep( f.Name, '\ |:', '_') ) );
        fprintf( '\n duration %.1f s (%.2f min)', toc - t, ( toc - t ) / 60 )
    end
end