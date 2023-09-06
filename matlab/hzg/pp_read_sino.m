%% Read sinogram %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    end
    sino_names_mat = NameCellToMat( sino_names );
    
    % Parameters
    filename = sprintf('%s%s', sino_path, sino_names_mat( round( par.im_shape_binned2 / 2 ), :));
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
    num_proj_used = num_proj_read;
    
    % Angles
    if  ~exist( 'angles', 'var' )
        if isempty( tomo.rot_angle_full_range )
            cprintf( 'Red', '\nEnter full angle of rotation (including one additional increment) or vector of angles, in radians: ' );
            tomo.rot_angle_full_range = input( '' );
        end
        if isscalar( tomo.rot_angle_full_range )
            angles = tomo.rot_angle_full_range * (0:num_proj_read - 1) / num_proj_read;
        else
            angles = tomo.rot_angle_full_range;
        end
        if length( angles ) ~= num_proj_read
            error( 'Number of angles (%u) entered not consistent with sinogram (%u) read.', numel( angles), num_proj_read )
        end
    end
    
    % Preallocation
    fprintf( '\nRead sinograms.')
    proj = zeros( par.im_shape_cropbin1, par.im_shape_binned2, num_proj_read, 'single');
    fprintf( ' Allocated bytes: %.2f GiB.', Bytes( proj, 3 ) )
    fprintf( ' shape : [num_pixel num_slice num_proj] = [%u %u %u]', par.im_shape_cropbin1, par.im_shape_binned2, num_proj_read )
    
    % Read sinogram    
    im_shape_cropbin1 = par.im_shape_cropbin1;    
    parfor nn = 1:size( proj, 2 )
        filename = sprintf('%s%s', sino_path, sino_names_mat(nn, :));
        sino = read_image( filename );
        sino = read_sino_trafo( sino );
        sino = sino(x0:x1,:);
        proj(:, nn, :) = reshape( sino, [im_shape_cropbin1, 1, num_proj_read] );
        %proj(:, nn, :) = permute( shiftdim( sino, -1 ) , [3 1 2] );
    end
    
    % Plot data
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
    fprintf( ' done in %.1f s (%.2f min)', toc - t, ( toc - t ) / 60 )
end