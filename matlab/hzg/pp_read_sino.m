%% Read sinogram script %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if par.read_sino
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
    else
        sino_names = {data_struct.name};
        par.im_shape_binned2 = numel( sino_names );
        im_shape_binned2 = par.im_shape_binned2;
    end
    sino_names_mat = NameCellToMat( sino_names );
    
    % Parameters
    filename = sprintf('%s%s', sino_path, sino_names_mat( round( par.im_shape_binned2 / 2 ),:));
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
    num_proj_read = size( sino, 1 );
    %num_proj_used = size( sino, 1 );
    num_proj_used = num_proj_read;
    
    
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
        
        % pixel size
        b = contains(c{1},'scan_pixelsize');
        if sum(b)
            pixelsize = c{2}(b);
            pixelsize = str2double(pixelsize{1});            
            fprintf('\n scan pixelsize: %f micron', pixelsize * 1000)
            c{1}(b) = [];
            c{2}(b) = [];
        end
        
        b = contains(c{1},'pixelsize');
        if sum(b)
            pixelsize = str2double(c{2}{b});
            fprintf('\n pixelsize: %f micron', pixelsize * 1000)
            eff_pixel_size = rawbin * pixelsize / 1000;
        end
        
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
    proj = zeros( par.im_shape_cropbin1, par.im_shape_binned2, num_proj_read, 'single');
    fprintf( '\n allocated memory: %.2f GiB.', Bytes( proj, 3 ) )
    fprintf( '\n shape: [num_pixel num_slice num_proj] = [%u %u %u]', par.im_shape_cropbin1, par.im_shape_binned2, num_proj_read )
    
    %% Read sinogram
    im_shape_cropbin1 = par.im_shape_cropbin1;
    im_shape_binned1 = im_shape_cropbin1;
    parfor nn = 1:size( proj, 2 )
        filename = sprintf('%s%s', sino_path, sino_names_mat(nn,:));
        sino = read_image( filename );
        sino = read_sino_trafo( sino );
        sino = sino(x0:x1,:);
        proj(:, nn,:) = reshape( sino, [im_shape_cropbin1, 1, num_proj_read] );
        %proj(:, nn,:) = permute( shiftdim( sino, -1 ) , [3 1 2] );
    end
    fprintf( '\n read in %.1f s (%.2f min)', toc - t, ( toc - t ) / 60 )
    
    %% Plot data
    if par.visual_output
        f = figure( 'Name', 'projection and sinogram', 'WindowState', par.window_state );
        
        subplot(1,2,1)
        imsc1( proj(:,:,1))
        title(sprintf('intensity: proj(:,:,1)'))
        colorbar
        axis equal tight
        
        subplot(1,2,2)
        nn = round( size( proj, 2) / 2);
        sino = squeeze( proj(:,nn,:) );
        sino = FilterOutlier( sino );
        imsc1( sino )
        title(sprintf('sinogram: proj(:,%u,:)', nn))
        colorbar
        axis equal tight
        
        drawnow
        
        CheckAndMakePath( write.fig_path )
        saveas( f, sprintf( '%s%s.png', write.fig_path, regexprep( f.Name, '\ |:', '_') ) );
    end
end