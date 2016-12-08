% P05 reconstruction pipeline. Preprocessing and tomographic
% reconstruction.
%
% USAGE
% Set parameters in PARAMETERS section and run script. If you are lucky you
% only have to adjust the 'scan_path' (or change to that folder and use
% 'pwd', see below).
%
% Written by Julian Moosmann.
% First version: 2016-09-28. Last modifcation: 2016-12-08

%clear all

%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%scan_path = pwd;
scan_path = ...
    '/asap3/petra3/gpfs/p05/2016/data/11001978/raw/mah_28_15R_top';
    '/asap3/petra3/gpfs/p05/2015/data/11001102/raw/hzg_wzb_mgag_14';    
    '/asap3/petra3/gpfs/p05/2016/data/11001978/raw/mah_28_15R_bottom';    
    '/asap3/petra3/gpfs/p05/2016/data/11001978/raw/mah_24_50L_top_load';
    '/asap3/petra3/gpfs/p05/2016/data/11001978/raw/mah_16_57R_load';
    '/asap3/petra3/gpfs/p05/2016/data/11001978/raw/mah_15_57R';
    '/asap3/petra3/gpfs/p05/2016/data/11001978/raw/mah_10_13R_bottom';
    '/asap3/petra3/gpfs/p05/2016/data/11001994/raw/szeb_74_13';
    '/asap3/petra3/gpfs/p05/2016/commissioning/c20160803_001_pc_test/raw/phase_1000';
    '/asap3/petra3/gpfs/p05/2016/commissioning/c20160803_001_pc_test/raw/phase_1400';
    '/asap3/petra3/gpfs/p05/2016/commissioning/c20160913_000_synload/raw/mg5gd_21_3w';
    '/asap3/petra3/gpfs/p05/2016/commissioning/c20160803_001_pc_test/raw/we43_phase_030';
    '/asap3/petra3/gpfs/p05/2016/data/11001464/raw/pnl_16_petrosia_c';
    '/asap3/petra3/gpfs/p05/2016/commissioning/c20160920_000_diana/raw/Mg-10Gd39_1w';
    '/asap3/petra3/gpfs/p05/2016/commissioning/c20161024_000_xeno/raw/xeno_01_b';
read_proj = 0; % Read flatfield-corrected images from disc
read_proj_folder = []; % subfolder of 'flat_corrected' containing projections
%proj_stride = 3; % stride of projection images to read
proj_range = 3; % range of found projections to be used. if empty: all, if scalar: stride
ref_range = [2:2:17]; % range of found flats to be used. if empty: all, if scalar: stride
bin = 4; % bin size: if 2 do 2 x 2 binning, if 1 do nothing
poolsize = 28; % number of workers in parallel pool to be used
gpu_ind = 1; % GPU Device to use: gpuDevice(gpu_ind)
gpu_thresh = 0.8; % Percentage of maximally used to available GPU memory
num_proj_hard = []; % hard coded number of projection. required if projections are missing
correct_beam_shake = 1;%  correlate flat fields and projection to correct beam shaking
correct_beam_shake_max_shift = 0; % if 0: use the best match (i.e. the one which is closest to zero), if > 0 all flats which are shifted less than correct_beam_shake_max_shift are used
flat_corr_area1 = [1 floor(100/bin)]; % correlation area: proper index range or relative/absolute position of [first pix, last pix]
flat_corr_area2 = [0.25 0.75]; %correlation area: proper index range or relative/absolute position of [first pix, last pix]
ring_filter = 1; % ring artifact filter
ring_filter_median_width = 11;
phase_retrieval = 0;
phase_retrieval_method = 'tie';
phase_retrieval_reg_par = 3;
phase_retrieval_bin_filt = 0.1;
energy = 30; % in keV
sample_detector_distance = 1.0; % in m
eff_pixel_size = 0.35e-6; % effective pixel size: (detector pixel size) / magnification
eff_pixel_size_binned = bin * eff_pixel_size; %1.3056e-5
do_tomo = 1; % reconstruct volume
vol_shape = []; % shape of the volume to be reconstructed, either in absolute number of voxels or in relative number w.r.t. the default volume which is given by the detector width and height
vol_size = []; % if empty, unit voxel size is assumed
rot_angle_full = []; % in radians: empty ([]), full angle of rotation, or array of angles. if empty full rotation angles is determined automatically
rot_angle_offset = pi; % global rotation of reconstructed volume
rot_axis_offset = []; % if empty use automatic computation
rot_axis_pos = []; % if empty use automatic computation. either offset or pos has to be empty. can't use both
rot_corr_area1 = [0.25 0.75]; % ROI to correlate projections at angles 0 & pi
rot_corr_area2 = [0.25 0.75]; % ROI to correlate projections at angles 0 & pi
rot_axis_tilt = 0 * -0.1 / 180 * pi; % camera tilt w.r.t rotation axis
fbp_filter_type = 'Ram-Lak';
fbp_filter_dim = 1; % dimension of sinogram to apply filte. should be the one of the detector lines orthogonal to the rotation axis
butterworth_filter = 1; % use butterworth filter in addition to FBP filter
butterworth_order = 1;
butterworth_cutoff_frequ = 0.5;
astra_pixel_size = 1; % size of a detector pixel: if different from one 'vol_size' needs to be ajusted 
link_data = 1; % ASTRA data objects become references to Matlab arrays.
take_neg_log = []; % logarithm for attenuation contrast
out_path = ''; % absolute path were output data will be stored. overwrites the write_to_scratch flage. if empty uses the beamtime directory and either 'processed' or 'scratch_cc'
write_proj = 0; % save preprocessed projections
write_reco = 1; % save reconstructed slices
write_to_scratch = 1; % write to 'scratch_cc' instead of 'processed'
parfolder = 'test'; % parent folder to 'reco' and 'flat_corrected'
parfolder_flatcor = ''; % parent folder to 'flat_corrected'
parfolder_reco = ''; % parent folder to 'reco'
verbose = 1; % print information to standard output
visualOutput = 0; % show images and plots during reconstruction

%% TODO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO: use log files instead of hard coded number of projection when
% images are missing
% TODO: merge subbranch to master
% TODO: automatic determination of rot center (entropy type)
% TODO: manual interactive finding of rotation center
% TODO: make pixel filtering thresholds parameters
% TODO: vertical ROI reco
% TODO: padding options for FBP filter
% TODO: normalize proj with beam current
% TODO: check offsets in projection correlation for rotation axis determination
% TODO: output file format option
% TODO: excentric rotation axis 
% TODO: save and read sinograms
% TODO: read number of projection from log file
% TODO: set photometric tag for tif files w/o one, turn on respective warning
% TODO: read image ROI and check if it's faster at all
% TODO: optional output format: 8-bit, 16-bit. Currently 32 bit tiff.
% TODO: check attenutation values of reconstructed slice
% TODO: ref stride
% TODO: adopt for missing images w.r.t. to img or tif
% TODO: save phase maps

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
warning( 'off', 'MATLAB:imagesci:rtifc:missingPhotometricTag');
%warning( 'off','all')
cdscandir = cd(scan_path);
PrintVerbose(verbose, '\nStart P05 reconstruction pipeline of scan_name: ')
NameCellToMat = @(name_cell) reshape(cell2mat(name_cell), [numel(name_cell{1}), numel(name_cell)])';
imsc1 = @(im) imsc( flipud( im' ) );
astra_clear % if reco was aborted, ASTRA memory is not cleared
if ~isempty( rot_axis_offset ) && ~isempty( rot_axis_pos )
    error('Either of rot_axis_offset (%f) or rot_axis_pos (%f) must be empty.', rot_axis_offset, rot_axis_pos)
end

%% Input folder
while scan_path(end) == filesep    
    scan_path(end) = [];
end
[raw_path, scan_name] = fileparts(scan_path);
scan_path = [scan_path, filesep];
[beamtime_dir, raw_folder] = fileparts(raw_path);
[~, beamtime_id] = fileparts(beamtime_dir);
if ~strcmp(raw_folder, 'raw')
    error('Name of folder is not raw: %s', raw_folder)
end
PrintVerbose(verbose, '%s', scan_name)
PrintVerbose(verbose, '\n raw_path:%s', scan_path)

%% Output folder
out_folder = 'processed'; 
if write_to_scratch(1)
    out_folder = 'scratch_cc';
end
if isempty( out_path )
    out_path = [beamtime_dir, filesep, out_folder, filesep, scan_name];
else
    out_path = [out_path, filesep, scan_name];
end
if ~isempty(parfolder)
    out_path = [out_path, filesep, parfolder];
end
% flat_corrected dir
if isempty( parfolder_flatcor )
    flatcor_path = [out_path, filesep, 'flat_corrected', filesep];
else
    flatcor_path = [out_path, filesep, 'flat_corrected', filesep, parfolder_flatcor, filesep];
end
PrintVerbose(verbose & write_proj, '\n flatcor: %s', flatcor_path)
% reco dir
if isempty( parfolder_reco )
    reco_dir = [out_path, filesep, 'reco', filesep];
else
    reco_dir = [out_path, filesep, 'reco', filesep, parfolder_reco, filesep];
end
PrintVerbose(verbose, '\n reco_path : %s', reco_dir)

%% File names

% Projection file names
proj_names = FilenameCell( [scan_path, '*.img'] );
if isempty( proj_names )
    proj_names =  FilenameCell( [scan_path, 'proj_*.tif'] );
end
proj_nums = CellString2Vec( proj_names );
num_proj_found = numel(proj_names);

% Ref file names
ref_names = FilenameCell( [scan_path, '*.ref'] );
if isempty( ref_names )
    ref_names = FilenameCell( [scan_path, 'ref_*.tif'] );
end
ref_nums = CellString2Vec( ref_names );
num_ref_found = numel(ref_names);
if isempty( ref_range )
    ref_range = 1;
end
if numel( ref_range ) == 1
    ref_range = 1:ref_range:num_ref_found;
end
num_ref_used = numel( ref_range );
PrintVerbose(verbose, '\n number of refs found : %g', num_ref_found)
PrintVerbose(verbose, '\n number of refs used : %g', num_ref_used)
PrintVerbose(verbose, '\n reference range used : %g:%g:%g%', ref_range(1), ref_range(2) - ref_range(1), ref_range(end))

% Dark file names
dark_names = FilenameCell( [scan_path, '*.dar'] );
if isempty( dark_names )
    dark_names = FilenameCell( [scan_path, 'dark_*.tif'] );
end
dark_nums = CellString2Vec( dark_names );
num_dark = numel(dark_names);
PrintVerbose(verbose, '\n number of darks found : %g', num_dark)

% Projection to read
if isempty( proj_range )
    proj_range = 1;
end
if numel( proj_range ) == 1
    proj_range = 1:proj_range:num_proj_found;
end
%proj_stride = max( 1, proj_stride);
%proj_range = 1:proj_stride:num_proj_found;
num_proj_used = numel( proj_range );
PrintVerbose(verbose, '\n number of projections found : %g', num_proj_found)
PrintVerbose(verbose, '\n number of projections used : %g', num_proj_used)
PrintVerbose(verbose, '\n projection range used : first:stride:last =  %g:%g:%g', proj_range(1), proj_range(2) - proj_range(1), proj_range(end))

%% Image shape and ROI
filename = sprintf('%s%s', scan_path, dark_names{1});
im = read_image( filename );
raw_im_shape = size( im );
raw_im_shape_binned = floor( raw_im_shape / bin );
raw_im_shape_binned1 = raw_im_shape_binned(1);
raw_im_shape_binned2 = raw_im_shape_binned(2);    
PrintVerbose(verbose, '\n raw image shape : %g  %g', raw_im_shape)
PrintVerbose(verbose, '\n raw image shape binned : %g  %g', raw_im_shape_binned)

%% Start parallel CPU pool
t = toc;
PrintVerbose( poolsize > 1, '\n Use parallel CPU pool of %u workers (or more). ', poolsize)
OpenParpool(poolsize);
PrintVerbose( poolsize > 1, ' Elapsed time: %.1f s', toc-t)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read flat corrected projection
if read_proj(1)
    data_struct = dir( [flatcor_path, 'proj*.*'] );
    if isempty( data_struct )
        
        fprintf('\n No flat corrected projections found! Switch to standard pre-processing.')
        read_proj = 0;
    else
        proj_names = {data_struct.name};
        num_proj_read = numel(proj_names);
        if num_proj_used > num_proj_read
            fprintf('\n Less projections available (%g) than demanded (%g)! Switch to standard pre-processing.', num_proj_read, num_proj_used )
            read_proj = 0;
        end
    end
    % File names
    t = toc;
    PrintVerbose(verbose, '\n Read flat corrected projections.')
    if num_proj_read ~= num_proj_used
        fprintf('\n CAUTION: Number of flat corrected projections read (%g) differs from number of projections to be processed (%g)!\n', num_proj_read, num_proj_used)
    end
    
    % Preallocation
    proj = zeros( raw_im_shape_binned(1), raw_im_shape_binned(2), num_proj_read, 'single');
    proj_names_mat = NameCellToMat( proj_names );
    
    % Read projections
    parfor nn = 1:num_proj_read
        filename = sprintf('%s%s', flatcor_path, proj_names_mat(nn, :));
        %proj(:, :, nn) = imread( filename, 'tif' )';
        proj(:, :, nn) = read_image( filename )';
    end
    PrintVerbose(verbose, ' Elapsed time: %.1f s (%.2f min)', toc - t, ( toc - t ) / 60 )
    
%% Read raw data
elseif ~read_proj
    %% Dark field
    t = toc;
    PrintVerbose(verbose, '\nProcessing dark fields.')
    dark = zeros( [raw_im_shape_binned, num_dark], 'single');
    parfor nn = 1:num_dark
        filename = sprintf('%s%s', scan_path, dark_names{nn});
        dark(:, :, nn) = Binning( FilterPixel( read_image( filename ), [0.02 0.02]), bin) / bin^2;
    end
    dark_min = min( dark(:) );
    dark_max = max( dark(:) );
    dark = squeeze( median(dark, 3) );
    %dark = FilterPixel( squeeze( median(dark, 3) ), [0.02 0.02]);
    if sum(dark <= 0)
        fprintf('\n CAUTION: dark field contains zeros')
    end
    PrintVerbose(verbose, ' Elapsed time: %.1f s', toc-t)
    if visualOutput(1)
        h1 = figure('Name', 'mean dark field, flat field, projections');
        subplot(2,2,1)
        imsc( dark' );
        axis equal tight;% square tight; pause(0.05);         
        title(sprintf('dark field'))
        drawnow
    end

    %% Flat field
    t = toc;
    PrintVerbose(verbose, '\nProcessing flat fields.')        
    
    % Preallocation
    flat = zeros( [raw_im_shape_binned, num_ref_used], 'single');
    ref_names_mat = NameCellToMat(ref_names);
    
    % Parallel loop
    parfor nn = 1:num_ref_used
        filename = sprintf('%s%s', scan_path, ref_names_mat(nn, :));        
        flat(:, :, nn) = Binning( FilterPixel( read_image( filename ), [0.01 0.005]), bin) / bin^2;
    end
    flat_min = min( flat(:) );
    flat_max = max( flat(:) );    
    if sum( flat(:) < 1 )
        fprintf('\n CAUTION: mean flat field contains zeros')
    end
    
    % Dark field correction
    flat = bsxfun( @minus, flat, dark );
    PrintVerbose(verbose, ' Elapsed time: %.1f s', toc-t)
    if visualOutput(1)
        figure(h1)
        subplot(2,2,2)
        imsc1( flat(:,:,1) );
        axis equal tight;% square
        title(sprintf('flat field #1'));
        drawnow;
    end

    %% Projections
    t = toc;
    PrintVerbose(verbose, '\nRead and filter raws.')
    % Preallocation
    proj = zeros( raw_im_shape_binned(1), raw_im_shape_binned(2), num_proj_used, 'single');    
    img_names_mat = NameCellToMat( proj_names );
    img_names_mat = img_names_mat( proj_range, :);
    % Display first raw images
    if visualOutput(1)  
        figure(h1)
        filename = sprintf('%s%s', scan_path, img_names_mat(1, :));
        raw1 = Binning( FilterPixel( read_image( filename ), [0.02 0.01]), bin) / bin^2;
        subplot(2,2,3)       
        imsc1( raw1 + dark );
        axis equal tight;% square
        title(sprintf('raw projection #1'))
        drawnow
    end
    % Read raw projections    
    parfor nn = 1:num_proj_used
        filename = sprintf('%s%s', scan_path, img_names_mat(nn, :));                
        proj(:, :, nn) = Binning( FilterPixel( read_image( filename ), [0.02 0.01]), bin) / bin^2;
    end    
    raw_min = min( proj(:) );
    raw_max = max( proj(:) );
    % Dark field correction
    proj = bsxfun( @minus, proj, dark);
    % Flat field correction without correlation
    if ~correct_beam_shake
        flat_m = median(flat, 3);
        proj = bsxfun( @times, proj, flat_m);
    end    
    PrintVerbose(verbose, ' Elapsed time: %.1f s (%.2f min)', toc - t, ( toc - t ) / 60 )    

    %% Correlate shifted flat fields
    if correct_beam_shake    
        PrintVerbose(verbose, '\nCorrect beam shake.')
        % Correlation ROI
        flat_corr_area1 = IndexParameterToRange(flat_corr_area1, raw_im_shape_binned(1));
        flat_corr_area2 = IndexParameterToRange(flat_corr_area2, raw_im_shape_binned(2));
        flatroi = flat(flat_corr_area1, flat_corr_area2, :);
        projroi = proj(flat_corr_area1, flat_corr_area2, :);
        xshift = zeros( num_proj_used, num_ref_used);
        yshift = zeros( num_proj_used, num_ref_used);
        % Compute shift for each pair projection/flat field
        for ff = 1:num_ref_used
            flat_ff = flatroi(:,:,ff);
            parfor pp = 1:num_proj_used
                out = ImageCorrelation(projroi(:,:,pp), flat_ff, 0,0);
                xshift(pp, ff) = out.Xshift;
                yshift(pp, ff) = out.Yshift; % relevant shift
            end    
        end       
        if visualOutput(1)
            h2 = figure('Name', 'Corrlation of raw data and flat fields');
            [~, yshift_min_pos] =  min ( abs( yshift), [], 2);
            [~, xshift_min_pos] =  min ( abs( xshift), [], 2);
            m = 1; n = 2;
            subplot(m,n,1);
            plot( arrayfun(@(x) (yshift(x,yshift_min_pos(x))), 1:num_proj_used) ,'.')            
            axis equal tight square
            title(sprintf('minimal shift: vertical'))
            subplot(m,n,2);
            plot( arrayfun(@(x) (xshift(x,xshift_min_pos(x))), 1:num_proj_used) ,'.')            
            axis equal tight square
            title(sprintf('minimal shift: horizontal (unused)'))            
            drawnow
        end
        
        % use best match
        if correct_beam_shake_max_shift == 0            
            [~, pos] = min( abs( yshift ), [], 2 );
            parfor nn = 1:num_proj_used
                %[~, pos] = min( abs( yshift(nn, : ) ) );                
                proj(:, :, nn) = proj(:, :, nn) ./ flat(:, :, pos(nn));
            end
        % use all flats which are shifted less pixels than correct_beam_shake_max_shift
        elseif correct_beam_shake_max_shift > 0
            flat_count = zeros(1, num_proj_used);
            parfor nn = 1:num_proj_used
                vec = 1:num_ref_used;
                flat_ind = vec( abs( yshift(nn, :) ) < correct_beam_shake_max_shift );
                if isempty( flat_ind )
                    [~, flat_ind] = min( abs( yshift(nn, :) ) );
                end
                flat_count(nn) = numel(flat_ind);    
                flat_selected_mean = FilterPixel( squeeze( mean( flat(:, :, flat_ind), 3) ), [0.005 0.0025]);           
                proj(:, :, nn) = proj(:, :, nn) ./ flat_selected_mean;
            end
        else
            error('Value of maximum shift (%g) is not >= 0', correct_beam_shake_max_shift)
        end
        PrintVerbose(verbose, ' Elapsed time: %.1f s (%.2f min)', toc - t, ( toc - t ) / 60 )
        if correct_beam_shake_max_shift > 0
            PrintVerbose(verbose, '\n number of flats used per projection: [mean, min, max] = [%g, %g, %g]', mean( flat_count ), min( flat_count ), max( flat_count) )
        end
    end
    PrintVerbose(verbose, '\n sinogram size = [%g, %g, %g]', size( proj ) )
    if visualOutput(1)
        figure(h1)
        subplot(2,2,4)        
        imsc1( proj(:,:,1));
        axis equal tight
        title(sprintf('flat-&-dark corrected projection #1'))
        drawnow
    end

    %% Ring artifact filter
    if ring_filter
        t = toc;
        PrintVerbose(verbose, '\nFilter ring artifacts.')
        sino = squeeze(proj(round(raw_im_shape_binned1/2),:,:));
        proj_mean = mean( proj, 3);
        proj_mean_med = medfilt2( proj_mean, [ring_filter_median_width, 1], 'symmetric' );
        mask = proj_mean_med ./ proj_mean;      
        proj = bsxfun( @times, proj, mask);            
        PrintVerbose(verbose, ' Elapsed time: %.1f s (%.2f min)', toc-t, (toc-t)/60)
        if visualOutput(1)          
            h3 = figure('Name', 'Sinogram and ring filter');
            subplot(2,2,1)
            imsc( sino )            
            axis equal tight
            title(sprintf('sinogram unfiltered'))
            subplot(2,2,2)
            imsc( squeeze(proj(round(raw_im_shape_binned1/2),:,:)) )
            axis equal tight
            title(sprintf('sinogram filtered'))            
            subplot(2,2,3)
            imsc1( proj_mean )
            axis equal tight
            title(sprintf('mean projection'))
            subplot(2,2,4)
            imsc1( mask )
            axis equal tight
            title(sprintf('mask for normalization'))            
            drawnow        
        end
    end
    
    proj_min = min( proj(:) );
    proj_max = max( proj(:) );
    
    %% Write corrected projections
    if write_proj(1)
        t = toc;   
        PrintVerbose(verbose, '\nSave flat-corrected projections.')
        CheckAndMakePath( flatcor_path )
        % Projections
        parfor nn = 1:num_proj_used     
            filename = sprintf('%sproj_%06u.tif', flatcor_path, nn );
            write32bitTIFfromSingle(filename, proj(:, :, nn)' );
        end
        PrintVerbose(verbose, ' Elapsed time: %g .0f (%.2f min)', toc-t, (toc-t)/60)
    end 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  %% Rotation axis position and Tomgraphic parameters %%%%%%%%%%%%%%%%%%%
if do_tomo    
    t = toc;
    rot_corr_area1 = IndexParameterToRange( rot_corr_area1, raw_im_shape_binned1 );
    rot_corr_area2 = IndexParameterToRange( rot_corr_area2, raw_im_shape_binned2 );
    % if numel(rot_corr_area1) == 2
    %     rot_corr_area1 = floor( rot_corr_area1(1) * raw_im_shape_binned(1) ):ceil( rot_corr_area1(2) * raw_im_shape_binned(1) );
    % end
    % if numel(rot_corr_area2) == 2
    %     rot_corr_area2 = floor( rot_corr_area2(1) * raw_im_shape_binned(2) ):ceil( rot_corr_area2(2) * raw_im_shape_binned(2) );
    % end

    % Automatic determination of full rotation angle if rot_angle_full is empty
    if isempty( rot_angle_full )        
        im1 = proj( rot_corr_area1, rot_corr_area2, 1);
        im2 = proj( rot_corr_area1, rot_corr_area2, num_proj_used);
        [~, cm1] = ImageCorrelation( im1, im2); % large if full angle is  2 * pi
        [~, cm2] = ImageCorrelation( im1, flipud(im2)); % large if full angle is pi
        if max( cm1(:) ) > max( cm2(:) )
            rot_angle_full = 2 * pi;
        elseif max( cm1(:) ) < max( cm2(:) )
            rot_angle_full = pi;
        else
            error('Automatic determination of full angle of rotation via correlation of first and last (flipped) projection failed.')
        end
    end
    %if numel( angles ) == 1
    PrintVerbose(verbose, '\n full rotation angule: %g', rot_angle_full)
    angles = rot_angle_full * (0:num_proj_used - 1) / (num_proj_used - 1);
% use above for img files and below for tiff files which accounts for missing images 
%     if isempty( num_proj_hard )
%         angles = rot_angle_full * proj_nums / num_proj_used;
%     else
%         angles = rot_angle_full * proj_nums / num_proj_hard;
%     end
%     else
%         PrintVerbose(verbose, '\n  angular range:', max( angles ) - min( angles ))
%     end
    if numel( angles ) ~= num_proj_used
        error('Number of elements in array of angles (%g) unequal number of projections read (%g)', numel( angles ), num_proj_used)
    end
    % retrieve index at angles 0 and pi
    [val1, ind1] = min( abs( angles ));
    [val2, ind2] = min( abs( angles - pi ));
    % Correlate images
    im1 = proj( :, rot_corr_area2, ind1);
    im2 = flipud( proj( :, rot_corr_area2, ind2) );    
    PrintVerbose(verbose, '\n image correlation: #%g at index %g and angle %g pi and #%g at index %g and angle %g pi', proj_nums(ind1), ind1, val1 / pi, proj_nums(ind2), ind2, (val2 + pi) / pi )
    out = ImageCorrelation( im1, im2, 0, 0);
    PrintVerbose(verbose, '\n relative shift: %g, %g', out.Xshift, out.Yshift)
    rot_axis_offset_calc = round( out.Xshift / 2, 1);     
    rot_axis_pos_calc = raw_im_shape_binned1 / 2 + rot_axis_offset_calc;        
    PrintVerbose(verbose, '\n calulated rotation axis offset: %g, %g', rot_axis_offset_calc)
    % use calculated offset if both offset and position are empty
    if isempty( rot_axis_offset ) && isempty( rot_axis_pos )
        rot_axis_offset = rot_axis_offset_calc;
    end
    if isempty( rot_axis_pos )
        rot_axis_pos = raw_im_shape_binned1 / 2 + rot_axis_offset;
    end
    if isempty( rot_axis_offset )
        rot_axis_offset = rot_axis_pos - raw_im_shape_binned1 / 2;
    end
    im1c = RotAxisSymmetricCropping( im1, rot_axis_pos, 1);
    im2c = flipud(RotAxisSymmetricCropping( flipud(im2), rot_axis_pos, 1));
    if visualOutput(1)
        h4 = figure('Name','Projections at 0 and pi cropped symmetrically to rotation center');
        subplot(2,2,1)
        imsc1( im1c )
        axis equal tight
        title(sprintf('proj at 0'))
        subplot(2,2,2)
        imsc1( im2c )
        axis equal tight
        title(sprintf('proj at pi'))
        subplot(2,2,3)
        imsc1( abs(im1c - im2c) )
        axis equal tight
        title(sprintf('difference'))
        drawnow 
        nimplay(cat(3, im1c',im2c'))
    end
    
    %% GPU slab sizes
    if isempty( vol_shape )
        % default volume given by the detector width and height
        vol_shape = [raw_im_shape_binned1, raw_im_shape_binned1, raw_im_shape_binned2];
    else            
        if vol_shape(1) <=  1
            vol_shape(1) = round( 1 + vol_shape(1) * (raw_im_shape_binned1 - 1) );
        end
        if vol_shape(2) <=  1
            vol_shape(2) = round( 1 + vol_shape(2) * (raw_im_shape_binned1 - 1) );
        end
        if vol_shape(3) <=  1
            vol_shape(3) = round( 1 + vol_shape(3) * (raw_im_shape_binned2 - 1) );
        end
    end
    if isempty( vol_size )
        vol_size = [-vol_shape(1)/2, vol_shape(1)/2, -vol_shape(2)/2, vol_shape(2)/2, -vol_shape(3)/2, vol_shape(3)/2];        
    end
    num_proj_used =  size( proj, 3);
    % GPU memory required for reconstruction
    rec_mem = ( raw_im_shape_binned1 * num_proj_used + vol_shape(1) * vol_shape(2) ) * 4;
    % GPU memory limit
    gpu = gpuDevice( gpu_ind );
    num_sli = floor( gpu_thresh * gpu.AvailableMemory / rec_mem );
    % slabs required
    num_slabs = ceil( vol_shape(3) / num_sli );
    % Readjust number slices to make all of similar size
    num_sli = ceil(vol_shape(3) / num_slabs);
    subvol_shape = [vol_shape(1:2), num_sli];
    PrintVerbose(verbose, '\n rotation axis offset: %g', rot_axis_offset );
    PrintVerbose(verbose, '\n rotation axis position: %g', rot_axis_pos );    
    PrintVerbose(verbose, '\n shape of reconstruction volume: [%g, %g, %g]', vol_shape )
    PrintVerbose(verbose, '\n shape of reconstruction subvolume: [%g, %g, %g]', subvol_shape )
    PrintVerbose(verbose, '\n available gpu memory : %g MiB (%.2f%%) of %g MiB', gpu.AvailableMemory/10^6, 100*gpu.AvailableMemory/gpu.TotalMemory, gpu.TotalMemory/10^6)
    PrintVerbose(verbose, '\n memory required for reconstruction of a single slice: %g MiB (estimate)', rec_mem / 1024^2 )
    PrintVerbose(verbose, '\n memory of reconstructed volume: %g MiB', prod( vol_shape ) * 4 / 1024^2 )
    PrintVerbose(verbose, '\n maximum memory of subvolume: %g MiB', prod( subvol_shape ) * 4 / 1024^2 )
    PrintVerbose(verbose, '\n number of subvolume slabs: %g', num_slabs )
    PrintVerbose(verbose, '\n maximum number of slices per slab: %g', num_sli )    
end

%% Phase retrieval %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if phase_retrieval(1)
        t = toc;   
        if isempty( take_neg_log )
            take_neg_log = 0;
        end
        
        % Phase retieval filter
        edp = [energy, sample_detector_distance, eff_pixel_size_binned];
        [pf, pha_appendix] = PhaseFilter( phase_retrieval_method, raw_im_shape_binned, edp, phase_retrieval_reg_par, phase_retrieval_bin_filt, 'single');
          
        % reco phase dir
        if isempty( parfolder_reco )
            reco_phase_dir = [out_path, filesep, 'reco_phase_', pha_appendix, filesep];
        else
            reco_phase_dir = [out_path, filesep, 'reco_phase_', pha_appendix, filesep, parfolder_reco, filesep];
        end
        CheckAndMakePath( reco_phase_dir )
        PrintVerbose(verbose, '\n reco phase  : %s', reco_phase_dir)
        PrintVerbose(verbose, '\nPhase retrieval.')
                     
        % Projections
        parfor nn = 1:num_proj_used            
            proj(:,:,nn) = -real( ifft2( pf .* fft2( proj(:,:,nn) ) ) );            
        end
        PrintVerbose(verbose, ' Elapsed time: %g s (%.2f min)', toc-t, (toc-t)/60)
        PrintVerbose(verbose, '\n phase retrieval method : %s', phase_retrieval_method)
        PrintVerbose(verbose, '\n energy : %g keV', energy)
        PrintVerbose(verbose, '\n sample detector distance : %g m', sample_detector_distance)        
end

%% Tomographic reco
if do_tomo
    PrintVerbose(verbose, '\nTomographic reconstruction of %u slabs:', num_slabs)
    vol_min = NaN;
    vol_max = NaN;
    for nn = 1:num_slabs

        % Slice indices
        sli0 = 1 + ( nn - 1) * num_sli;
        sli1 = min( [nn * num_sli, vol_shape(3)] );

        % Filter sinogram
        PrintVerbose(verbose, '\n slab %g of %g: Filter sino: ', nn, num_slabs)                     
        filt = iradonDesignFilter('Ram-Lak', size(proj,1), 1);
        if butterworth_filter            
            [b, a] = butter(butterworth_order, butterworth_cutoff_frequ);
            bw = freqz(b, a, numel(filt) );
            filt = filt .* bw;
        end
        sino = real( ifft( bsxfun(@times, fft( NegLog(permute(proj(:,sli0:sli1,:), [1 3 2]), take_neg_log), [], 1), filt), [], 1, 'symmetric') );        
        PrintVerbose(verbose, 'done.')

        % Backprojection
        PrintVerbose(verbose, ' Backproject: ')    
        subvol_shape(3) = size( sino, 3);
        subvol_size = vol_size;
        subvol_size(5) = subvol_shape(3) / vol_shape(3) * vol_size(5);
        subvol_size(6) = subvol_shape(3) / vol_shape(3) * vol_size(6);
        vol = astra_parallel3D( sino, rot_angle_offset + angles, rot_axis_offset, subvol_shape, subvol_size, astra_pixel_size, link_data, rot_axis_tilt);
        PrintVerbose(verbose, 'done.')

        % Save subvolume
        if write_reco
            PrintVerbose(verbose, ' Save slices:')
            if phase_retrieval(1)
                save_dir = reco_phase_dir;
            else
                save_dir = reco_dir;
            end
            CheckAndMakePath(reco_dir)
            % Save reco path to reconstruction in file
            filename = [userpath, filesep, 'experiments/p05/pathtolastreco'];
            fid = fopen( filename , 'w' );
            fprintf( fid, '%s', reco_dir );   
            fclose( fid );
            % Save slices
            parfor ii = 1:size( vol, 3)
                filename = sprintf( '%sreco_%06u.tif', save_dir, sli0 + ii - 1);
                write32bitTIFfromSingle( filename, vol( :, :, ii) )
            end
            PrintVerbose(verbose, ' done.')
        end
        vol_min = min( min( vol(:) ), vol_min );
        vol_max = max( max( vol(:) ), vol_max );
    end
    PrintVerbose(verbose, ' Elapsed time: %.1f s (%.2f min)', toc-t, (toc-t)/60 )
    
    % Write log file
    if write_reco
         filename = sprintf( '%sreco.log', save_dir);
         fid = fopen(filename, 'w');
         fprintf(fid, 'scan_name : %s\n', scan_name);
         fprintf(fid, 'beamtime_id : beamtime_id\n');
         fprintf(fid, 'scan_path : %s\n', scan_path);
         fprintf(fid, 'reco_path : %s\n', save_dir);                  
         fprintf(fid, 'MATLAB notation. Index of first element: 1. Range: first:stride:last\n');
         fprintf(fid, 'MATLAB version: %s\n', version);
         fprintf(fid, 'platform: %s\n', computer);
         fprintf(fid, 'num_dark_found : %u\n', num_dark);
         fprintf(fid, 'num_ref_found : %u\n', num_ref_found);
         fprintf(fid, 'num_ref_used : %u\n', num_ref_used);
         fprintf(fid, 'ref_range : %u:%u:%u\n', ref_range(1), ref_range(2) - ref_range(1), ref_range(end) );         
         fprintf(fid, 'num_proj_found : %u\n', num_proj_found);
         fprintf(fid, 'num_proj_used : %u\n', num_proj_used);
         fprintf(fid, 'proj_range : %u\n', proj_range(1), proj_range(2) - proj_range(1), proj_range(end) );
         fprintf(fid, 'raw_image_shape : %u %u\n', raw_im_shape);
         fprintf(fid, 'raw_image_shape_binned : %u %u\n', raw_im_shape_binned);
         fprintf(fid, 'raw_image_binning : %u\n', bin);
         fprintf(fid, 'effective_pixel_size : %g\n', eff_pixel_size);
         fprintf(fid, 'effective_pixel_size_binned : %g\n', eff_pixel_size_binned);
         fprintf(fid, 'energy_eV : %g\n', energy);
         fprintf(fid, 'flat_field_correlation_area_1: %u:%u:%u\n', flat_corr_area1(1), flat_corr_area1(2) - flat_corr_area1(1), flat_corr_area1(end));
         fprintf(fid, 'flat_field_correlation_area_2: %u:%u:%u\n', flat_corr_area2(1), flat_corr_area2(2) - flat_corr_area2(1), flat_corr_area2(end));
         fprintf(fid, 'min_max_of_all_darks : %6g %6g\n', dark_min, dark_max);         
         fprintf(fid, 'min_max_of_all_flats : %6g %6g\n', flat_min, flat_max);
         fprintf(fid, 'min_max_of_all_raws :  %6g %6g\n', raw_min, raw_max);         
         fprintf(fid, 'min_max_of_all_flat_corr_projs : %g %g \n', proj_min, proj_max);
         % Phase retrieval
         fprintf(fid, 'use_phase_retrieval : %u\n', phase_retrieval);
         fprintf(fid, 'phase_retrieval_method : %s\n', phase_retrieval_method);
         fprintf(fid, 'phase_retrieval_regularisation_parameter : %f\n', phase_retrieval_reg_par);
         fprintf(fid, 'phase_retrieval_binary_filter_threshold : %f\n', phase_retrieval_bin_filt);
         % Rotation
         fprintf(fid, 'rotation_angle_full : %f\n', rot_angle_full);
         fprintf(fid, 'rotation_angle_offset : %f\n', rot_angle_offset);
         fprintf(fid, 'rotation_axis_offset_calculated : %f\n', rot_axis_offset_calc);
         fprintf(fid, 'rotation_axis_offset_used : %f\n', rot_axis_offset);
         fprintf(fid, 'rotation_axis_position_calculated : %f\n', rot_axis_pos_calc);
         fprintf(fid, 'rotation_axis_position_used : %f\n', rot_axis_pos);
         fprintf(fid, 'raw_image_binned_center : %f\n', raw_im_shape_binned1 / 2);         
         fprintf(fid, 'rotation_correlation_area_1 : %u:%u:%u\n', rot_corr_area1(1), rot_corr_area1(2) - rot_corr_area1(1), rot_corr_area1(end));
         fprintf(fid, 'rotation_correlation_area_2 : %u:%u:%u\n', rot_corr_area2(1), rot_corr_area2(2) - rot_corr_area2(1), rot_corr_area2(end));         
         % FBP
         fprintf(fid, 'use_ring_filter : %u\n', ring_filter);
         fprintf(fid, 'ring_filter_median_width : %u\n', ring_filter_median_width);
         fprintf(fid, 'fbp_filter : %s\n', fbp_filter_type);
         fprintf(fid, 'use_butterworth_filter : %u\n', butterworth_filter);
         fprintf(fid, 'butterworth_order : %u\n', butterworth_order);
         fprintf(fid, 'butterworth_cutoff_frequency : %f\n', butterworth_cutoff_frequ);
         fprintf(fid, 'astra_pixel_size : %f\n', astra_pixel_size);
         fprintf(fid, 'take_negative_logarithm : %u\n', take_neg_log);
         fprintf(fid, 'gpu_name : %s\n', gpu.Name);
         fprintf(fid, 'num_subvolume_slabs : %u\n', num_slabs);
         fprintf(fid, 'min_max_of_all_slices : %g %g\n', vol_min, vol_max);         
         %fprintf(fid, ' : \n', );
         fprintf(fid, 'seconds_elapsed_for_full_reco : %.1f\n', toc);
         fclose(fid);
         PrintVerbose(verbose, '\n log file : %s', filename)
    end

end

PrintVerbose(verbose, '\nFinished. Total time elapsed: %g s (%.2f min)\n\n', toc, toc / 60 );
% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
