%% Read flat corrected projections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if par.read_flatcor
    t = toc;
    read_flatcor_bin = par.read_flatcor_bin;
    read_flatcor_range = par.read_flatcor_range;
    par.raw_bin = read_flatcor_bin;
    
    if isempty( par.read_flatcor_path )
        par.read_flatcor_path = flatcor_path;
    end
    
    %     % Angles
    %     filename = [ write.reco_path filesep 'angles.mat' ];
    %     if ~exist( filename, 'file' )
    %         p = [fileparts( fileparts( par.read_flatcor_path ) ) filesep 'reco'];
    %         filename = [ p filesep 'angles.mat' ];
    %         if ~exist( filename, 'file')
    %             p = [ fileparts( par.read_flatcor_path )  filesep 'reco'];
    %             filename = [ p filesep 'angles.mat' ];
    %         end
    %         if ~exist( filename, 'file')
    %             error( 'angles.mat not found' )
    %         end
    %     end
    %     load( filename, 'angles' );
    
    if isempty(par.sample_detector_distance)
        d = dir([par.nexus_path filesep '*.h5']);
        nexuslog_name = [d.folder filesep d.name];
        h5log_group = h5info(nexuslog_name, '/entry/scan/setup/' );
        if sum(strcmp('o_ccd_dist',{h5log_group.Datasets.Name}))
            par.sample_detector_distance = h5read( nexuslog_name, '/entry/scan/setup/o_ccd_dist');
        elseif sum(strcmp('pos_o_ccd_dist',{h5log_group.Datasets.Name}))
            par.sample_detector_distance = h5read( nexuslog_name, '/entry/scan/setup/pos_o_ccd_dist');
        else
            warning('5h entrey ''o_ccd_dist'' not found. Setting distance to 1')
            par.sample_detector_distance = 1;
        end
        par.sample_detector_distance = par.sample_detector_distance / 1000;
        fprintf('\n sample dector distance: %.3f m', par.sample_detector_distance )
    end
    if isempty(par.energy)
        d = dir([par.nexus_path filesep '*.h5']);
        nexuslog_name = [d.folder filesep d.name];
        
        h5log_group = h5info(nexuslog_name, '/entry/scan/setup/' );
        if sum(strcmp('p07_energy',{h5log_group.Datasets.Name}))
            par.energy= h5read( nexuslog_name, '/entry/scan/setup/p07_energy');
        elseif sum(strcmp('pos_p05_energy',{h5log_group.Datasets.Name}))
            par.energy= h5read( nexuslog_name, '/entry/scan/setup/pos_p05_energy');
        else
            warning('h5 entry ''o_ccd_dist entry'' not found . Setting energy to 1')
            par.energy = 1;
        end
        fprintf('\n energy: %.1f keV', par.energy / 1e3 )
    end
    if ~exist('angles','var') || isempty(angles)
        d = dir([par.nexus_path filesep '*.h5']);
        nexuslog_name = [d.folder filesep d.name];
        s_rot.value = h5read( nexuslog_name, '/entry/scan/data/s_rot/value');
        par.num_dark = h5read( nexuslog_name, '/entry/scan/n_dark');
        [~, stimg_key, ~, ~] = pp_stimg_petra({nexuslog_name},par);
        angles = s_rot.value( ~boolean( stimg_key.scan.value(par.num_dark+1:end) ) ) * pi / 180;
    end
    
    % File names
    data_struct = dir( [par.read_flatcor_path filesep '*.tif'] );
    if isempty( data_struct )
        fprintf('\n No flat corrected projections found! Switch to standard pre-processing.')
        par.read_flatcor = 0;
        %     else
        %         proj_names = {data_struct.name};
        %
        %         % Projection range to read
        %         if isempty( par.proj_range )
        %             par.proj_range = 1;
        %         end
        %         if numel( par.proj_range ) == 1
        %             par.proj_range = 1:par.proj_range:num_proj_found;
        %         end
        %         num_proj_used = numel( par.proj_range );
        %         %proj_nums = CellString2Vec( proj_names(par.proj_range), imtype_str_flag);
        %         fprintf( '\n projections found : %g', num_proj_found)
        %         fprintf( '\n projections used : %g', num_proj_used)
        %         fprintf( '\n projection range used : first:stride:last =  %g:%g:%g', par.proj_range(1), par.proj_range(2) - par.proj_range(1), par.proj_range(end))
    end
    flatcor_names = {data_struct.name};
    num_proj_found = numel(flatcor_names);
    if isempty(read_flatcor_range)
        read_flatcor_range = 1;
    end
    if numel(read_flatcor_range) == 1
        read_flatcor_range = 1:read_flatcor_range:num_proj_found;
    end
    flatcor_names_used = flatcor_names(read_flatcor_range);
    angles = angles(read_flatcor_range);
    num_proj_used = numel(flatcor_names_used);
    fprintf( '\n projections found : %g', num_proj_found)
    fprintf( '\n projections used : %g', num_proj_used)
    if ~isempty(tomo.rot_angle_full_range)
        if isscalar( tomo.rot_angle_full_range )
            angles = tomo.rot_angle_full_range * (0:num_proj_used - 1) / num_proj_used;
        else
            angles = tomo.rot_angle_full_range;
        end
    end
    
    fprintf( '\nReading flat corrected projections: ' )
    % Preallocation
    filename = sprintf('%s/%s', par.read_flatcor_path, flatcor_names_used{1});
    read_flatcor_trafo = par.read_flatcor_trafo;
    im = read_image( filename );
    im = Binning(im,read_flatcor_bin);
    im = read_flatcor_trafo( im );
    proj = zeros( [ size( im ), num_proj_used ], 'single' );
    fprintf( '\n allocated bytes: %.2f GiB.', Bytes( proj, 3 ) )
    flatcor_names_mat = NameCellToMat( flatcor_names_used );
    % Read flat corrected projections
    read_flatcor_path = par.read_flatcor_path;
    parfor nn = 1:num_proj_used
        filename = sprintf('%s/%s', read_flatcor_path, flatcor_names_mat(nn, :));
        im = read_image(filename);
        im = Binning(im,read_flatcor_bin);
        im = read_flatcor_trafo(im);
        proj(:, :, nn) = im;
    end
    fprintf( '\n done in %.1f s (%.2f min)', toc - t, ( toc - t ) / 60 )
end
