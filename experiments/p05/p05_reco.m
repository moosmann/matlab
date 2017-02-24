% P05 reconstruction pipeline. Preprocessing, filtering, phase retrieval,
% and tomographic reconstruction.
%
% USAGE
% Set parameters in PARAMETERS / SETTINGS section below and run script.
%
% To start script: type 'F5' when focus is in the Editor windows, or click
% 'Run' in the Editor tab, or type 'p05_reco_loop' and Enter in the Command
% Window. To loop over data or parameter use 'p05_reco_loop'.
%
% Written by Julian Moosmann. First version: 2016-09-28. Last modifcation:
% 2017-02-22

%% PARAMETERS / SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%scan_path = pwd; % Enter folder of data set under the raw directory and run script
scan_path = ...
    '/asap3/petra3/gpfs/p05/2016/commissioning/c20161024_000_xeno/raw/xeno_01_b';
    '/asap3/petra3/gpfs/p05/2016/commissioning/c20160803_001_pc_test/raw/phase_600';
    '/asap3/petra3/gpfs/p05/2016/commissioning/c20160803_001_pc_test/raw/phase_1000'; %rot_axis_offset=7.5    
    '/asap3/petra3/gpfs/p05/2016/commissioning/c20160803_001_pc_test/raw/phase_1400'; %rot_axis_offset=19.5        
    '/asap3/petra3/gpfs/p05/2016/commissioning/c20160913_000_synload/raw/mg5gd_21_3w';
    '/asap3/petra3/gpfs/p05/2016/data/11001994/raw/szeb_41';
    '/asap3/petra3/gpfs/p05/2016/commissioning/c20160803_001_pc_test/raw/we43_phase_030';    
    '/asap3/petra3/gpfs/p05/2016/data/11001978/raw/mah_32_15R_top_occd800_withoutpaper'; % no rings in FS
    '/asap3/petra3/gpfs/p05/2016/data/11001978/raw/mah_30_15R_top_occd125_withoutpaper'; % no rings in FS
    '/asap3/petra3/gpfs/p05/2016/data/11001994/raw/szeb_23_00'; % rot_axis_offset 5.75;
    '/asap3/petra3/gpfs/p05/2016/data/11001994/raw/szeb_23_01'; % Nothobranchius furzeri
    '/asap3/petra3/gpfs/p05/2016/data/11001994/raw/szeb_23_00'; % Nothobranchius furzeri; bewegung
    '/asap3/petra3/gpfs/p05/2016/data/11001994/raw/szeb_13_00'; % no conspicuous movement artifacts, but cell shape are unclear and nuclei not visible
    '/asap3/petra3/gpfs/p05/2016/data/11001994/raw/szeb_13_09'; % dead    
    '/asap3/petra3/gpfs/p05/2016/data/11001994/raw/szeb_07_00'; % strong movement
    '/asap3/petra3/gpfs/p05/2016/data/11001994/raw/szeb_13_00';
    '/asap3/petra3/gpfs/p05/2015/data/11001102/raw/hzg_wzb_mgag_14';        
    '/asap3/petra3/gpfs/p05/2016/data/11001978/raw/mah_36_1R_top';    
    '/asap3/petra3/gpfs/p05/2016/data/11001978/raw/mah_28_15R_top';    
    '/asap3/petra3/gpfs/p05/2015/data/11001102/raw/hzg_wzb_mgag_38';
    '/asap3/petra3/gpfs/p05/2015/data/11001102/raw/hzg_wzb_mgag_02';           
    '/asap3/petra3/gpfs/p05/2016/data/11001978/raw/mah_28_15R_bottom';    
    '/asap3/petra3/gpfs/p05/2016/data/11001978/raw/mah_24_50L_top_load';
    '/asap3/petra3/gpfs/p05/2016/data/11001978/raw/mah_16_57R_load';
    '/asap3/petra3/gpfs/p05/2016/data/11001978/raw/mah_15_57R';
    '/asap3/petra3/gpfs/p05/2016/data/11001978/raw/mah_10_13R_bottom';
    '/asap3/petra3/gpfs/p05/2016/data/11001464/raw/pnl_16_petrosia_c';
    '/asap3/petra3/gpfs/p05/2016/commissioning/c20160920_000_diana/raw/Mg-10Gd39_1w';    
% Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
read_flatcor = 0; % Read flatfield-corrected images from disc. Skips preprocessing
read_flatcor_path = '/asap3/petra3/gpfs/p05/2016/data/11001978/scratch_cc/c20160803_001_pc_test/phase_1000/flat_corrected'; % subfolder of 'flat_corrected' containing projections
% Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out_path = '/gpfs/petra3/scratch/moosmanj'; % absolute path were output data will be stored. !!overwrites the write_to_scratch flag. if empty uses the beamtime directory and either 'processed' or 'scratch_cc'
%out_path = '/asap3/petra3/gpfs/p05/2016/data/11001978/scratch_cc/c20160803_001_pc_test';
write_to_scratch = 0; % write to 'scratch_cc' instead of 'processed'
write_flatcor = 1; % save preprocessed flat corrected projections
write_phase_map = 1; % save phase maps (if phase retrieval is not 0)
write_sino = 1; % save sinograms (after preprocessing & before FBP filtering and phase retrieval)
write_sino_phase = 1; % save sinograms of phase maps
write_reco = 1; % save reconstructed slices (if do_tomo=1)
parfolder = ''; % parent folder of 'reco', 'sino', and 'flat_corrected'
subfolder_flatcor = ''; % subfolder of 'flat_corrected'
subfolder_phase_map = ''; % subfolder of 'phase_map'
subfolder_sino = ''; % subfolder of 'sino'
subfolder_reco = '';%sprintf('fbpFilt%s_ringFilt%uMedWid%u_bwFilt%ubwCutoff%u_phasePad%u_freqCutoff%2.0f_fbpPad%u', fbp_filter_type, ring_filter, ring_filter_median_width, butterworth_filter, 100*butterworth_cutoff_frequ, phase_padding, fpb_freq_cutoff*100, fbp_filter_padding); % parent folder to 'reco'
% Preprocessing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bin = 2; % projection binning factor: 1, 2, or 4
do_stitching = 1; % stitch projection of 2 pi scan at rotation axis position
stitching_method = 'sine'; % 'step': no interpolation, 'linear','sine': linear interpolation of overlap region. !!! adjust: correlation area
proj_range = 1; % range of found projections to be used. if empty: all, if scalar: stride
ref_range = 1; % range of flat fields to be used: start:inc:end. if empty: all (equals 1). if scalar: stride
darkFiltPixHot = 0.01; % Hot pixel filter parameter for dark fields, for details see 'FilterPixel'
darkFiltPixDark = 0.005; % Dark pixel filter parameter for dark fields, for details see 'FilterPixel'
refFiltPixHot = 0.01; % Hot pixel filter parameter for flat fields, for details see 'FilterPixel'
refFiltPixDark = 0.005; % Dark pixel filter parameter for flat fields, for details see 'FilterPixel'
projFiltPixHot = 0.01; % Hot pixel filter parameter for projections, for details see 'FilterPixel'
projFiltPixDark = 0.005; % Dark pixel filter parameter for projections, for details see 'FilterPixel'
correct_beam_shake = 1;%  correlate flat fields and projection to correct beam shaking
correct_beam_shake_max_shift = 0.5; % if 0: use the best match (i.e. the one with the least shift), if > 0 all flats which are shifted less than correct_beam_shake_max_shift are used
max_num_flats = 9; % maximum number of flat fields used for average/median of flats
norm_by_ring_current = 1; % normalize flat fields and projections by ring current
flat_corr_area1 = [1 floor(100/bin)]; % correlation area: proper index range or relative/absolute position of [first pix, last pix]
flat_corr_area2 = [0.25 0.75]; %correlation area: proper index range or relative/absolute position of [first pix, last pix]
round_precision = 2; % precision when rounding pixel shifts
ring_filter = 1; % ring artifact filter
ring_filter_median_width = 11;
% Phase retrieval %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
do_phase_retrieval = 1;
phase_retrieval_method = 'tie';'qpcut'; %'qp' 'ctf' 'tie' 'qp2'
phase_retrieval_reg_par = 2.5; % regularization parameter
phase_retrieval_bin_filt = 0.15; % threshold for quasiparticle retrieval 'qp', 'qp2'
phase_retrieval_cutoff_frequ = 1 * pi; % in radian. frequency cutoff in Fourier space for 'qpcut' phase retrieval
phase_padding = 1; % padding of intensities before phase retrieval
energy = []; % in eV. if empty: read from log file
sample_detector_distance = []; % in m. if empty: read from log file
eff_pixel_size = []; % in m. if empty: read from log file. effective pixel size =  detector pixel size / magnification
% Tomography parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
do_tomo = 1; % run tomographic reconstruction of volume
vol_shape = []; % shape of the volume to be reconstructed, either in absolute number of voxels or in relative number w.r.t. the default volume which is given by the detector width and height
vol_size = []; % if empty, unit voxel size is assumed
rot_angle_full = []; % in radians: empty ([]), full angle of rotation, or array of angles. if empty full rotation angles is determined automatically to pi or 2 pi
rot_angle_offset = pi; % global rotation of reconstructed volume
rot_axis_offset = [];19.5;7.5;5.75; % if empty use automatic computation
rot_axis_pos = []; % if empty use automatic computation. either offset or pos has to be empty. can't use both
rot_corr_area1 = [0.75 1]; % ROI to correlate projections at angles 0 & pi. Use [0.75 1] or so for scans with an excentric rotation axis
rot_corr_area2 = [0.2 0.8]; % ROI to correlate projections at angles 0 & pi
rot_axis_tilt = 0 * -0.1 / 180 * pi; % camera tilt w.r.t rotation axis
fbp_filter_type = 'linear';'Ram-Lak';
fpb_freq_cutoff = 1; % Cut-off frequency in Fourier space of the above FBP filter
fbp_filter_padding = 1; % symmetric padding for consistent boundary conditions
butterworth_filter = 1; % use butterworth filter in addition to FBP filter
butterworth_order = 1;
butterworth_cutoff_frequ = 0.5;
astra_pixel_size = 1; % size of a detector pixel: if different from one 'vol_size' needs to be ajusted 
take_neg_log = []; % take negative logarithm. if empty, use 1 for attenuation contrast, 0 for phase contrast
% Interaction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
verbose = 1; % print information to standard output
visualOutput = 0; % show images and plots during reconstruction
interactive_determination_of_rot_axis = 1; % reconstruct slices with different offsets
interactive_determination_of_rot_axis_slice = 0.5; % slice number. if in [0,1]: relative, if in (1, N]: absolute
% Hardware / Software %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
poolsize = 0.9; % number of workers used in a parallel pool. if > 1: absolute number. if 0 < poolsize < 1: relative amount of all cores to be used
gpu_ind = 1; % GPU Device to use: gpuDevice(gpu_ind)
gpu_thresh = 0.8; % Percentage of maximally used to available GPU memory
image_corr_gpu = 0; % experimental, slower than CPU: uses GPU to correlate images
link_data = 1; % ASTRA data objects become references to Matlab arrays.

%% TODO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO: stitching: optimize and refactor, memory efficency, exponentian interpolatin method
% TODO: interactive loop over tomo slices for different phase retrieval parameter
% TODO: automatic determination of rot center: entropy type
% TODO: output file format option: 8-bit, 16-bit. Currently 32 bit tiff.
% TODO: normalize proj with beam current for KIT camera and missing images
% TODO: vertical ROI reco
% TODO: read image ROI and test for speed up
% TODO: additional padding schemes for FBP filter
% TODO: read sinograms option
% TODO: set photometric tag for tif files w/o one, turn on respective warning
% TODO: parallelize slice reconstruction using parpool and astra
% TODO: GPU phase retrieval: parfor-loop requires memory managment
% TODO: check attenutation values of reconstructed slice are consitent
% TODO: median filter width of ring filter dependence on binning
% TODO: check offsets in projection correlation for rotation axis determination
% TODO: delete files before writing data to a folder
% TODO: log file location for non-tomo processing
% TODO: inverse Gaussian filter for phase retrieval, VZC theorem
% TODO: check calculated rotation axis offset when horizontal ROI is used for image correlation

%% Notes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Padding and phase retrieval:
% Symmetric padding of intensity maps before phase retrieval cleary reduces
% artifacts in the retrieved phase maps which are due to inconistent
% boundary conditions. However, due to local tomography phase retrieval
% with padding can result in more artifacts in the recostructed volume
% (such as additonal bumps at the halo described below). After phase
% retrieval without padding the region close to the left
% boundary is blended in to the region close to the right boundary and vice
% versa (same holds for top and bottom). This effect fades out with
% increasing distance to the boundarys. Due to this blending effect 
% asymmetrical 'density' distributions, which are the cause of halo-like
% artifacts stretching outwards from the biggest possible circle within the
% reconstruction volume within a slice, are reduced which reduces the
% halo-like artifacts as well.

% FOV extension by excentric rotation axis:
% For absorpion-contrast data (more precisely when no phase retrieval
% is used), volumes can be reconstructed from a data set where an excentric
% position of the rotation axis is used without prior stitching of the
% projections by simply providing the correct rotation axis position and
% setting fbp_filter_padding to 1. For the automatic detection of the
% rotation axis to work, the area to correlate has to be adjusted i.e.
% 'rot_corr_area1' (not yet tested). When a phase retrieval is used, this
% approach does not work appropriately and gives rise to artifacts near the
% center of the reconstructed volume. This is due to the fact, that without
% stitching the phase is retrieved from a 'cropped' projection which
% results in inconsitently retrieved low frequencies (large scale
% variations) in the phase map. Using the 'linear' FBP filter instead of
% 'Ram-Lak' can maybe reduce these artifacts (not tested).

%% External call: parameters set by 'p05_reco_loop' %%%%%%%%%%%%%%%%%%%%%%%
if exist( 'external_parameter' ,'var')
    visualOutput = 0; 
    interactive_determination_of_rot_axis = 0;
    fields = fieldnames( external_parameter );
    for nn = 1:numel(fields)
        var_name = fields{nn};
        var_val = getfield(external_parameter, var_name );
        assignin('caller', var_name, var_val )
    end    
    clear external_parameter;
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
PrintVerbose(verbose, 'Start reconstruction of ')
warning( 'off', 'MATLAB:imagesci:rtifc:missingPhotometricTag');
cdscandir = cd(scan_path);
NameCellToMat = @(name_cell) reshape(cell2mat(name_cell), [numel(name_cell{1}), numel(name_cell)])';
imsc1 = @(im) imsc( flipud( im' ) );
astra_clear % if reco was aborted, ASTRA memory is not cleared
if ~isempty( rot_axis_offset ) && ~isempty( rot_axis_pos )
    error('rot_axis_offset (%f) and rot_axis_pos (%f) cannot be used simultaneously. One must be empty.', rot_axis_offset, rot_axis_pos)
end
if do_tomo || image_corr_gpu
    gpu = gpuDevice( gpu_ind );
    reset( gpu );
    gpu = gpuDevice( gpu_ind );
end

%% Input folder
while scan_path(end) == filesep    
    scan_path(end) = [];
end
[raw_path, scan_name] = fileparts(scan_path);
scan_path = [scan_path, filesep];
[beamtime_path, raw_folder] = fileparts(raw_path);
[~, beamtime_id] = fileparts(beamtime_path);
if ~strcmp(raw_folder, 'raw')
    error('Name of folder is not raw: %s', raw_folder)
end
PrintVerbose(verbose, '%s', scan_name)
PrintVerbose(verbose, '\n scan_path:%s', scan_path)
% Save scan path in file            
filename = [userpath, filesep, 'experiments/p05/pathtolastscan'];
fid = fopen( filename , 'w' );
fprintf( fid, '%s', scan_path );
fclose( fid );

%% Output folder
out_folder = 'processed'; 
if write_to_scratch(1)
    out_folder = 'scratch_cc';
end
if isempty( out_path )
    out_path = [beamtime_path, filesep, out_folder, filesep, scan_name];
else
    out_path = [out_path, filesep, scan_name];
end
if ~isempty(parfolder)
    out_path = [out_path, filesep, parfolder];
end

% Projections path
if isempty( subfolder_flatcor )
    flatcor_path = [out_path, filesep, 'flat_corrected', filesep];
else
    flatcor_path = [out_path, filesep, 'flat_corrected', filesep, subfolder_flatcor, filesep];
end
PrintVerbose(verbose & write_flatcor, '\n flatcor_path: %s', flatcor_path)

% Phase map path
if isempty( subfolder_phase_map )
       phase_map_path = [out_path, filesep, 'phase_map', filesep];
else
    phase_map_path = [out_path, filesep, 'phase_map', filesep, subfolder_phase_map, filesep];    
end
PrintVerbose(verbose & write_phase_map, '\n phase_map_path: %s', phase_map_path)

% Sinogram path
if isempty( subfolder_sino )
    sino_path = [out_path, filesep, 'sino', filesep];
    sino_phase_path = [out_path, filesep, 'sino_phase', filesep];
else
    sino_path = [out_path, filesep, 'sino', filesep, subfolder_sino, filesep]; 
    sino_phase_path = [out_path, filesep, 'sino_phase', filesep, subfolder_sino, filesep]; 
end
PrintVerbose(verbose & write_sino, '\n sino_path: %s', sino_path)
PrintVerbose(verbose & write_sino_phase, '\n sino_phase_path: %s', sino_phase_path)

% Reco path
if isempty( subfolder_reco )
    reco_path = [out_path, filesep, 'reco', filesep];
else
    reco_path = [out_path, filesep, 'reco', filesep, subfolder_reco, filesep];
end
PrintVerbose(verbose, '\n reco_path : %s', reco_path)

%% File names

% Projection file names
proj_names = FilenameCell( [scan_path, '*.img'] );
if isempty( proj_names )
    proj_names =  FilenameCell( [scan_path, 'proj_*.tif'] );
end
num_proj_found = numel(proj_names);

% Ref file names
ref_names = FilenameCell( [scan_path, '*.ref'] );
if isempty( ref_names )
    ref_names = FilenameCell( [scan_path, 'ref_*.tif'] );
end
num_ref_found = numel(ref_names);
if isempty( ref_range )
    ref_range = 1;
end
if numel( ref_range ) == 1
    ref_range = 1:ref_range:num_ref_found;
end
num_ref_used = numel( ref_range );
ref_names_mat = NameCellToMat( ref_names(ref_range) );
ref_nums = CellString2Vec( ref_names(ref_range) );
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

% Projection range to read
if isempty( proj_range )
    proj_range = 1;
end
if numel( proj_range ) == 1
    proj_range = 1:proj_range:num_proj_found;
end
num_proj_used = numel( proj_range );
proj_nums = CellString2Vec( proj_names(proj_range) );
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

%% P05 log-file
str = dir( sprintf( '%s*scan.log', scan_path) );
filename = sprintf( '%s/%s', str.folder, str.name);
[par, cur, cam] = p05_log( filename );
if isempty( energy )
    if isfield( par, 'Energy' )
        energy = par.Energy;
    elseif isfield( par, 'energy' )
        energy = par.energy;
    end              
end
if isempty( eff_pixel_size )
    if isfield( par, 'eff_pix_size' )
        eff_pixel_size = par.eff_pix_size * 1e-3 ;
    elseif isfield( par, 'eff_pix' )
        eff_pixel_size = par.eff_pix * 1e-3 ;
    elseif isfield( par, 'ccd_pixsize' ) && isfield( par, 'magn' )
        eff_pixel_size = par.ccd_pixsize / par.magn * 1e-3 ;
    end
end
eff_pixel_size_binned = bin * eff_pixel_size;
if isempty( sample_detector_distance )
    if isfield( par, 'camera_distance')
        sample_detector_distance = par.camera_distance / 1000;
    elseif isfield( par, 'o_ccd_dist')
        sample_detector_distance = par.o_ccd_dist / 1000;
    end
end

%% Start parallel CPU pool
t = toc;
if poolsize <= 1 && poolsize > 0
    poolsize = max( floor( poolsize * feature('numCores') ), 1 );
end
PrintVerbose( poolsize > 1, '\nStart parallel pool of %u workers. ', poolsize)
OpenParpool(poolsize);
PrintVerbose( poolsize > 1, ' Time elapsed: %.1f s', toc-t)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read flat corrected projection
if read_flatcor(1)
    t = toc;
    if isempty( read_flatcor_path )
        read_flatcor_path = flatcor_path;
    end
    % File names    
    data_struct = dir( [read_flatcor_path filesep 'proj*.*'] );
    if isempty( data_struct )        
        fprintf('\n No flat corrected projections found! Switch to standard pre-processing.')
        read_flatcor = 0;
    else
        proj_names = {data_struct.name};
        num_proj_read = numel(proj_names);
        if num_proj_used > num_proj_read
            fprintf('\n Less projections available (%g) than demanded (%g)! Switch to standard pre-processing.', num_proj_read, num_proj_used )
            read_flatcor = 0;
        end
    end    
    PrintVerbose(verbose, '\n Read flat corrected projections.')
    if num_proj_read ~= num_proj_used
        fprintf('\n WARNING: Number of flat corrected projections read (%g) differs from number of projections to be processed (%g)!\n', num_proj_read, num_proj_used)
    end
    
    % Preallocation
    proj = zeros( raw_im_shape_binned(1), raw_im_shape_binned(2), num_proj_read, 'single');
    proj_names_mat = NameCellToMat( proj_names );
    
    % Read flat corrected projections
    parfor nn = 1:num_proj_read
        filename = sprintf('%s%s', flatcor_path, proj_names_mat(nn, :));        
        proj(:, :, nn) = fliplr( read_image( filename ) );
    end
    PrintVerbose(verbose, ' Time elapsed: %.1f s (%.2f min)', toc - t, ( toc - t ) / 60 )
    
%% Read raw data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ~read_flatcor
    %% Dark field   
    t = toc;
    PrintVerbose(verbose, '\nProcessing dark fields.')
    dark = zeros( [raw_im_shape_binned, num_dark], 'single');
    parfor nn = 1:num_dark
        filename = sprintf('%s%s', scan_path, dark_names{nn});
        im = single( read_image( filename ) );
        % Remove large outliers. Assume Poisson distribtion at large lambda
        % is approximately a Gaussian distribution and set all value above
        % mean + 4 * std (99.994 of values lie within 4 std). Due to
        % outliers 4*std will contain much more values and is a good
        % estimate
        im_mean = mean( im(:) );
        im_std = std( im(:) );
        im( im > im_mean + 4*im_std) = im_mean;
        dark(:, :, nn) = Binning( FilterPixel( im, [darkFiltPixHot darkFiltPixDark]), bin) / bin^2;
    end
    dark_min = min( dark(:) );
    dark_max = max( dark(:) );
    dark = squeeze( median(dark, 3) );
    dark_med_min = min( dark(:) );
    dark_med_max = max( dark(:) );
    PrintVerbose(verbose, ' Time elapsed: %.1f s', toc-t)
    if visualOutput(1)
        h1 = figure('Name', 'mean dark field, flat field, projections');
        subplot(2,2,1)
        imsc1( dark );
        axis equal tight;% square tight; pause(0.05);         
        title(sprintf('dark field'))
        colorbar
        drawnow
    end

    %% Flat field
    t = toc;
    PrintVerbose(verbose, '\nProcessing flat fields.')        
    
    % Preallocation
    flat = zeros( [raw_im_shape_binned, num_ref_used], 'single');    
    
    % Parallel loop
    parfor nn = 1:num_ref_used
        filename = sprintf('%s%s', scan_path, ref_names_mat(nn, :));        
        flat(:, :, nn) = Binning( FilterPixel( read_image( filename ), [refFiltPixHot refFiltPixDark]), bin) / bin^2;        
        % Check for zeros        
        num_zeros =  sum( sum( flat(:,:,nn) < 1 ) );
        if num_zeros > 0
            fprintf('\n WARNING: values of %u pixels of flat field no %u are lower than 1.', num_zeros, nn)            
        end        
    end
    % min/max values before dark field subtraction and ring normalization
    flat_min = min( flat(:) );
    flat_max = max( flat(:) ); 
    
    % Dark field correction
    flat = bsxfun( @minus, flat, dark );
    
    % Ring current normalization
    if norm_by_ring_current(1)
        if isequal( ref_nums, [cur.ref(ref_range).ind] )
            scale_factor = 100 ./ shiftdim( [cur.ref(ref_range).val], -1 );
            flat = bsxfun( @times, flat, scale_factor );
            if visualOutput(1)
                hrc = figure('Name', 'Ring current normalization');
                subplot(1,2,1);
                plot( scale_factor(:) )
                axis equal tight square
                title(sprintf('flats'))
                drawnow
            end            
        else
            fprintf('\n WARNING: flat fields not normalized by ring current. Names read from dir() and log-file are inconsistent.')
        end
    end
    
    % Check for non-positive values
    parfor nn = 1:num_ref_used
        im = flat(:,:,nn);
        m = im < 1;
        im(m) = 1;
        flat(:,:,nn) = im;
    end
    flat_min2 = min( flat(:) );
    flat_max2 = max( flat(:) ); 
    
    num_zeros =  sum( flat(:) < 1 );
    if num_zeros > 0
        fprintf('\n WARNING: flat field contains %u zeros', num_zeros)        
    end
        
    PrintVerbose(verbose, ' Time elapsed: %.1f s', toc-t)
    if visualOutput(1)
        figure(h1)
        subplot(2,2,2)
        imsc1( flat(:,:,1) )
        axis equal tight% square
        title(sprintf('flat field #1'))
        drawnow
        colorbar
    end

    %% Projections
    t = toc;
    PrintVerbose(verbose, '\nRead and filter raws.')
    % Preallocation
    proj = zeros( raw_im_shape_binned(1), raw_im_shape_binned(2), num_proj_used, 'single');    
    img_names_mat = NameCellToMat( proj_names(proj_range) );
    
    % Display first raw images
    if visualOutput(1)  
        figure(h1)
        filename = sprintf('%s%s', scan_path, img_names_mat(1, :));
        raw1 = Binning( FilterPixel( read_image( filename ), [projFiltPixHot, projFiltPixDark]), bin) / bin^2;
        subplot(2,2,3)       
        imsc1( raw1 )
        axis equal tight% square
        title(sprintf('raw projection #1'))
        drawnow
        colorbar
    end
    
    % Read raw projections
    parfor nn = 1:num_proj_used
        filename = sprintf('%s%s', scan_path, img_names_mat(nn, :));                
        proj(:, :, nn) = Binning( FilterPixel( read_image( filename ), [projFiltPixHot, projFiltPixDark]), bin) / bin^2;
    end    
    raw_min = min( proj(:) );
    raw_max = max( proj(:) );
    
    % Dark field correction
    proj = bsxfun( @minus, proj, dark);
    
    % Ring current normalization
    if norm_by_ring_current(1)
        if isequal( proj_nums, [cur.proj(proj_range).ind] )
            scale_factor = 100 ./ shiftdim( [cur.proj(proj_range).val], -1 );
            proj = bsxfun( @times, proj, scale_factor );            
            if visualOutput(1)
                figure(hrc)
                subplot(1,2,2);
                plot( scale_factor(:) )
                axis equal tight square
                title(sprintf('projections'))
                drawnow
            end            
        else
            fprintf('\n WARNING: projections not normalized by ring current. Names read from dir() and log-file are not consistent.')
        end
    end    
    
    % Check for non-positive values
    parfor nn = 1:num_proj_used
        im = proj(:,:,nn);
        m = im < 1;
        if sum( m(:) ) > 0            
            im(m) = 1;
            proj(:,:,nn) = im;
        end
    end
    raw_min2 = min( proj(:) );
    raw_max2 = max( proj(:) );
    
    %% Flat field correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Flat field correction without correlation
    if ~correct_beam_shake
        flat_m = median(flat, 3);
        proj = bsxfun( @times, proj, flat_m);
    end    
    PrintVerbose(verbose, ' Time elapsed: %.1f s (%.2f min)', toc - t, ( toc - t ) / 60 )    

    % Correlate shifted flat fields
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
            
            if image_corr_gpu(1) % GPU, no parloop
                for pp = 1:num_proj_used
                    out = ImageCorrelation(projroi(:,:,pp), flat_ff, 0, 0, 1);
                    xshift(pp, ff) = round( out.Xshift, round_precision );
                    yshift(pp, ff) = round( out.Yshift, round_precision) ; % relevant shift
                end                
            else % CPU, parloop
                parfor pp = 1:num_proj_used
                    out = ImageCorrelation(projroi(:,:,pp), flat_ff, 0, 0, 0);
                    xshift(pp, ff) = round( out.Xshift, round_precision );
                    yshift(pp, ff) = round( out.Yshift, round_precision) ; % relevant shift
                end
            end            
        end       
        if visualOutput(1)
            h2 = figure('Name', 'Correlation of raw data and flat fields');
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
                proj(:, :, nn) = proj(:, :, nn) ./ flat(:, :, pos(nn));
            end
        % use all flats which are shifted less pixels than correct_beam_shake_max_shift
        elseif correct_beam_shake_max_shift > 0
            flat_count = zeros(1, num_proj_used);
            parfor nn = 1:num_proj_used
                vec = 1:num_ref_used;
                flat_ind = vec( abs( yshift(nn, :) ) < correct_beam_shake_max_shift );
                if numel( flat_ind ) > max_num_flats
                    flat_ind( max_num_flats + 1:end ) = [];
                end
                if isempty( flat_ind )
                    [~, flat_ind] = min( abs( yshift(nn, :) ) );
                end
                flat_count(nn) = numel(flat_ind);                    
                proj(:, :, nn) = proj(:, :, nn) ./ squeeze( mean( flat(:, :, flat_ind), 3) );
            end
        else
            error('Value of maximum shift (%g) is not >= 0', correct_beam_shake_max_shift)
        end
        PrintVerbose(verbose, ' Time elapsed: %.1f s (%.2f min)', toc - t, ( toc - t ) / 60 )
        if correct_beam_shake_max_shift > 0
            PrintVerbose(verbose, '\n number of flats used per projection: [mean, min, max] = [%g, %g, %g]', mean( flat_count ), min( flat_count ), max( flat_count) )
        end
    end
    PrintVerbose(verbose, '\n sinogram size = [%g, %g, %g]', size( proj ) )
    if visualOutput(1)
        figure(h1)
        subplot(2,2,4)        
        imsc1( proj(:,:,1))
        axis equal tight
        title(sprintf('flat-&-dark corrected projection #1'))
        colorbar
        drawnow        
    end

    %% Ring artifact filter
    if ring_filter(1)
        t = toc;
        PrintVerbose(verbose, '\nFilter ring artifacts.')        
        proj_mean = mean( proj, 3);
        proj_mean_med = medfilt2( proj_mean, [ring_filter_median_width, 1], 'symmetric' );
        mask = proj_mean_med ./ proj_mean;      
        proj = bsxfun( @times, proj, mask);            
        PrintVerbose(verbose, ' Time elapsed: %.1f s (%.2f min)', toc-t, (toc-t)/60)
        PrintVerbose( verbose, '\n ring filter mask min/max: %f, %f', min( mask(:) ), max( mask(:) ) )
        if visualOutput(1)
            sino = squeeze(proj(round(raw_im_shape_binned1/2),:,:));
            h3 = figure('Name', 'Sinogram and ring filter');
            subplot(2,2,1)
            imsc( sino' )            
            axis equal tight
            title(sprintf('sinogram unfiltered'))
            colorbar('Location', 'southoutside')
            subplot(2,2,2)
            imsc( squeeze(proj(round(raw_im_shape_binned1/2),:,:))' )
            axis equal tight
            title(sprintf('sinogram filtered'))
            colorbar('Location', 'southoutside')
            subplot(2,2,3)
            imsc1( proj_mean )
            axis equal tight
            title(sprintf('mean projection'))
            colorbar
            subplot(2,2,4)
            imsc1( mask )
            axis equal tight
            title(sprintf('mask for normalization'))
            colorbar
            drawnow        
        end
    end    
    proj_min = min( proj(:) );
    proj_max = max( proj(:) );
    
    %% Write corrected projections
    if write_flatcor(1)
        t = toc;   
        PrintVerbose(verbose, '\nSave flat-corrected projections.')
        CheckAndMakePath( flatcor_path )
        % Projections
        parfor nn = 1:num_proj_used     
            filename = sprintf('%sproj_%06u.tif', flatcor_path, nn );
            write32bitTIFfromSingle(filename, rot90( proj(:, :, nn) ) );
        end
        PrintVerbose(verbose, ' Time elapsed: %g .0f (%.2f min)', toc-t, (toc-t)/60)
    end 
end

%% Rotation axis position and tomgraphic reconstruction parameters %%%%%%%%
if do_tomo(1)
    t = toc;    
    rot_corr_area1 = IndexParameterToRange( rot_corr_area1, raw_im_shape_binned1 );
    rot_corr_area2 = IndexParameterToRange( rot_corr_area2, raw_im_shape_binned2 );

    % Full rotation angle
    if isempty( rot_angle_full )
        if isfield( par, 'rotation')
            % From log file
            rot_angle_full = par.rotation / 180 * pi;
        elseif exist('cur', 'var') && isfield(cur, 'proj') && isfield( cur.proj, 'angle')
            % from beam current log
            rot_angle_full = (cur.proj(end).angle - cur.proj(1).angle) * pi /180; % KIT: , EHD: ok
        else
            % Guess from correlation of projections
            im1 = proj( rot_corr_area1, rot_corr_area2, 1);
            im2 = proj( rot_corr_area1, rot_corr_area2, num_proj_used);
            [~, cm1] = ImageCorrelation( im1, im2, 0, 0, image_corr_gpu); % large if full angle is  2 * pi
            [~, cm2] = ImageCorrelation( im1, flipud(im2), 0, 0, image_corr_gpu); % large if full angle is pi
            if max( cm1(:) ) > max( cm2(:) )
                rot_angle_full = 2 * pi;
            elseif max( cm1(:) ) < max( cm2(:) )
                rot_angle_full = pi;
            else
                error('Determination of full angle of rotation from correlation of first and last projection not successful.')
            end
        end
    end    
    PrintVerbose(verbose, '\n full rotation angle: %g * pi', rot_angle_full / pi)
    if exist('cur', 'var') && isfield(cur, 'proj') && isfield( cur.proj, 'angle')
        % for KIT cam this includes missing angles
        angles = [cur.proj.angle] / 180 * pi; 
        if strcmpi(cam, 'kit')
            % drop angles where projections are missing
            angles = angles(1 + proj_nums);
        else
            angles = angles(proj_range); 
        end
    else
        if isfield( par, 'projections' )
            num_proj = double( par.projections );
        elseif isfield( par, 'n_angles' )
            num_proj = double( par.n_angles );
        end        
        switch lower( cam )
            case 'ehd'
                angles = rot_angle_full * (0:num_proj - 1) / (num_proj - 1); % EHD: ok
            case 'kit'
                angles = rot_angle_full * (0:num_proj - 1) / num_proj; % KIT: ok if par.projections exist
        end        
    end    
    if numel( angles ) ~= num_proj_used
        error('Number of elements in array of angles (%g) unequal number of projections read (%g)', numel( angles ), num_proj_used)
    end
    % retrieve index at angles 0 and pi
    [val1, ind1] = min( abs( angles ));
    [val2, ind2] = min( abs( angles - pi ));
    % Correlate images
%     im1 = proj( :, rot_corr_area2, ind1);
%     im2 = flipud( proj( :, rot_corr_area2, ind2) );    
    im1 = proj( rot_corr_area1, rot_corr_area2, ind1);
    im2 = flipud( proj( rot_corr_area1, rot_corr_area2, ind2) );    
    PrintVerbose(verbose, '\n image correlation: #%g at index %g and angle %g pi and #%g at index %g and angle %g pi', proj_nums(ind1), ind1, val1 / pi, proj_nums(ind2), ind2, (val2 + pi) / pi )
    out = ImageCorrelation( im1, im2, 0, 0, image_corr_gpu);
    xshift = round( out.Xshift, round_precision) + rot_corr_area1(1) - 1;
    yshift = round( out.Yshift, round_precision) + rot_corr_area2(1) - 1;    
    PrintVerbose(verbose, '\n relative shift: %g, %g', xshift, yshift)
    rot_axis_offset_calc = xshift / 2;
    rot_axis_pos_calc = raw_im_shape_binned1 / 2 + rot_axis_offset_calc;        
    PrintVerbose(verbose, '\n calulated rotation axis offset: %f', rot_axis_offset_calc)
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
        colorbar
        subplot(2,2,2)
        imsc1( im2c )
        axis equal tight
        title(sprintf('proj at pi'))
        colorbar
        subplot(2,2,3)
        imsc1( abs(im1c - im2c) )
        axis equal tight
        title(sprintf('difference'))
        colorbar
        drawnow 
        nimplay(cat(3, im1c',im2c'))
    end

    %% Determine rotation axis position
    if interactive_determination_of_rot_axis(1)
        fprintf( '\n\nENTER INTERACTIVE MODE' )
        % default range is centered at the given or calculated offset
        offset = rot_axis_offset + (-4:0.5:4);
        fprintf( '\n calculated rotation axis offset & position : %.2f, %.2f', rot_axis_offset_calc, rot_axis_pos_calc)
        fprintf( '\n current rotation axis offset & position    : %.2f, %.2f', rot_axis_offset, rot_axis_pos)
        
        % Loop over offsets
        while ~isscalar( offset )
            pause(0.2)
            if interactive_determination_of_rot_axis_slice > 1
                slice = interactive_determination_of_rot_axis_slice;
            elseif interactive_determination_of_rot_axis_slice <= 1 && interactive_determination_of_rot_axis_slice >= 0
                slice = round((raw_im_shape_binned2 - 1) * interactive_determination_of_rot_axis_slice + 1 );
                %slice = floor(raw_im_shape_binned2 / 2);
            end
            [vr, m1, m2, m3, m4, m5, com] = find_rot_axis(proj, angles, offset, slice);
            fprintf( '\n sinogram center of mass offset : %g', com - raw_im_shape_binned1)
            
            % Print image number, rotation axis values, and different metrics
            fprintf( '\n image_no    offset         m1         m2         m3         m4         m5')
            for nn = 1:numel(offset)
                fprintf( '\n %8u %9.2f %10.3g %10.3g %10.3g  %10.3g %10.3g', nn, offset(nn), m1(nn), m2(nn), m3(nn), m4(nn), m5(nn) )
            end        
            nimplay(vr)            
            offset = input( '\nTYPE ROTATION AXIS OFFSET OR ENTER TO CONTINUE SCRIPT OR TYPE OFFSET RANGE AS [...] TO RE-RUN RECONSTRUCTIONS: ');
            if isempty( offset )
                fprintf( ' new and old rotation axis offset : %.2f', rot_axis_offset)
                break
            elseif isscalar( offset )                             
                fprintf( ' old rotation axis offset : %.2f', rot_axis_offset)
                rot_axis_offset = offset;
                fprintf( '\n new rotation axis offset : %.2f', rot_axis_offset)
                break
            else                
                fprintf( '\n new rotation axis positions and metrics :')
            end
        end
        rot_axis_pos = raw_im_shape_binned1 / 2 + rot_axis_offset;
    end
end

%% Stitch projections
if do_stitching(1)
    t = toc;
    PrintVerbose(verbose, '\nStitch projections:')
    if rot_angle_full < 1.9 * pi
        error( 'full angle of rotation smaller than 2 pi: %g pi', rot_angle_full/pi)
    end
    % number of unstitched projections
    num_proj = size( proj, 3);
    % last projection within [0,pi)
    [~, num_proj_sti] = min( abs(angles - pi));
    % number of stitched projections
    num_proj_sti = num_proj_sti - 1;
    % index range of projections to be stitched
    xl = 1:round(rot_axis_pos);
    xr = 1:xl(end)-1;
    im_shape_sti1 = numel( xl ) + numel( xr );
    % Preallocation
    proj_sti = zeros( im_shape_sti1 , raw_im_shape_binned2, num_proj_sti, 'single');    
    for nn = 1:num_proj_sti
        nn2 = mod(num_proj_sti + nn - 1, num_proj) + 1;        
        im = zeros( im_shape_sti1, raw_im_shape_binned2);
        switch lower( stitching_method )
            case 'step'
                im = cat(1, proj(xl,:,nn), flipud( proj(xr,:,nn2) ) );
            case {'linear', 'sine'}
                % overlap region
                overlap = 2 * rot_axis_pos - raw_im_shape_binned1 : raw_im_shape_binned1;
                % overlap ramp
                x = (0:1/(numel(overlap)-1):1);
                % 1D weight
                w = ones(raw_im_shape_binned1, 1);
                switch lower( stitching_method )
                    case 'linear'
                        w(overlap) = 1 - x;
                    case 'sine'
                        w(overlap) = 0.5 * cos(pi*x) + 0.5;
                end                
                % weighted projections
                iml = bsxfun(@times, proj(:,:,nn), w);
                imr = flipud( bsxfun(@times, proj(:,:,nn2), w ) );
                % stitched projection
                im(1:raw_im_shape_binned1,:) = iml;
                im(end - raw_im_shape_binned1 + 1:end,:) = im(end - raw_im_shape_binned1 + 1:end,:) + imr;                
        end                
        proj_sti(:,:,nn) = im;
    end
    proj = proj_sti;
    clear proj_sti;
    angles = angles(1:num_proj_sti);
    PrintVerbose(verbose, ' Time elapsed: %g .0f (%.2f min)', toc-t, (toc-t)/60)
end

%% Save sinograms
if write_sino(1)  
    t = toc;
    PrintVerbose(verbose, '\nSave sinogram:')    
    CheckAndMakePath(sino_path)    
    % Save slices
    parfor nn = 1:raw_im_shape_binned2
        filename = sprintf( '%ssino_%06u.tif', sino_path, nn);
        write32bitTIFfromSingle( filename, squeeze( proj( :, nn, :) )' )
    end
    PrintVerbose(verbose, ' done.')
    PrintVerbose(verbose, ' Time elapsed: %g s (%.2f min)', toc-t, (toc-t)/60)   
end

%% Phase retrieval %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
edp = [energy, sample_detector_distance, eff_pixel_size_binned];
if do_phase_retrieval(1)
    t = toc;
    PrintVerbose(verbose, '\n energy : %g eV', energy)
    PrintVerbose(verbose, '\n sample detector distance : %g m', sample_detector_distance)
    
    if isempty( take_neg_log )
        take_neg_log = 0;
    end
    
    % Phase retrieval filter
    im_shape = [size(proj,1) , size(proj,2)];
    im_shape_pad = (1 + phase_padding) * im_shape;
    [phase_filter, pha_appendix] = PhaseFilter( phase_retrieval_method, im_shape_pad, edp, phase_retrieval_reg_par, phase_retrieval_bin_filt, phase_retrieval_cutoff_frequ, 'single');
           
    % reco phase dir
    if isempty( subfolder_reco )
        reco_phase_path = [out_path, filesep, 'reco_phase', filesep, pha_appendix, filesep];
    else
        reco_phase_path = [out_path, filesep, 'reco_phase', filesep, pha_appendix, filesep, subfolder_reco, filesep];
    end
    CheckAndMakePath( reco_phase_path )
    PrintVerbose(verbose, '\n reco_phase_path : %s', reco_phase_path)
    PrintVerbose(verbose, '\n phase retrieval method : %s', phase_retrieval_method)
    PrintVerbose(verbose, '\nPhase retrieval.')
    
    % Retrieval
    parfor nn = 1:size(proj, 3)
        % combined GPU and parfor usage requires memory management
        if phase_padding
            im = padarray( proj(:,:,nn), im_shape, 'symmetric', 'post' );
            %im = padarray( gpuArray( proj(:,:,nn) ), raw_im_shape_binned, 'post', 'symmetric' );
            im = -real( ifft2( phase_filter .* fft2( im ) ) );
            proj(:,:,nn) = im(1:im_shape(1), 1:im_shape(2));
            %proj(:,:,nn) = gather( im(1:raw_im_shape_binned1, 1:raw_im_shape_binned2) );
        else
            proj(:,:,nn) = -real( ifft2( phase_filter .* fft2( proj(:,:,nn) ) ) );
            %proj(:,:,nn) = gather( -real( ifft2( phase_filter .* fft2( gpuArray( proj(:,:,nn) ) ) ) ) );
        end
    end
    PrintVerbose(verbose, ' Time elapsed: %g s (%.2f min)', toc-t, (toc-t)/60)
    
    % Save phase maps
    if write_phase_map(1)
        t = toc;
        PrintVerbose(verbose, '\nSave phase maps:')
        phase_map_path = [phase_map_path, pha_appendix, filesep];
        CheckAndMakePath( phase_map_path );
        parfor nn = 1:size( proj, 3)
            filename = sprintf( '%sphase_%06u.tif', phase_map_path, nn);
            write32bitTIFfromSingle( filename, squeeze( rot90( proj( :, :, nn) ) ) )
        end
        PrintVerbose(verbose, ' Time elapsed: %g s (%.2f min)', toc-t, (toc-t)/60)
    end
end

%% Save sinograms of phase maps
if write_sino(1)
    t = toc;
    PrintVerbose(verbose, '\nSave phase map sinogram:')    
    CheckAndMakePath(sino_phase_path)    
    % Save slices
    parfor nn = 1:raw_im_shape_binned2
        filename = sprintf( '%ssino_%06u.tif', sino_phase_path, nn);
        write32bitTIFfromSingle( filename, squeeze( proj( :, nn, :) )' )
    end
    PrintVerbose(verbose, ' done.')
    PrintVerbose(verbose, ' Time elapsed: %g s (%.2f min)', toc-t, (toc-t)/60)   
end

%% Tomographic reco
if do_tomo(1)
    
    if do_stitching(1)
        rot_axis_offset_reco = 0;
    else
        rot_axis_offset_reco = rot_axis_offset;
    end
    
    %% GPU slab sizes
    if isempty( vol_shape )
        % default volume given by the detector width and height
        vol_shape1 = size( proj, 1);
        vol_shape = [vol_shape1, vol_shape1, raw_im_shape_binned2];
    else            
        if vol_shape(1) <=  1
            vol_shape(1) = round( 1 + vol_shape(1) * (vol_shape1 - 1) );
        end
        if vol_shape(2) <=  1
            vol_shape(2) = round( 1 + vol_shape(2) * (vol_shape1 - 1) );
        end
        if vol_shape(3) <=  1
            vol_shape(3) = round( 1 + vol_shape(3) * (raw_im_shape_binned2 - 1) );
        end
    end
    if isempty( vol_size )
        vol_size = [-vol_shape(1)/2, vol_shape(1)/2, -vol_shape(2)/2, vol_shape(2)/2, -vol_shape(3)/2, vol_shape(3)/2];        
    end    
    % GPU memory required for reconstruction of single slice
    rec_mem = ( vol_shape1 * size( proj, 3) + vol_shape(1) * vol_shape(2) ) * 4;
    % GPU memory limit    
    num_sli = floor( gpu_thresh * gpu.AvailableMemory / rec_mem );
    % slabs required
    num_slabs = ceil( vol_shape(3) / num_sli );
    % Readjust number slices to make all of similar size
    num_sli = ceil(vol_shape(3) / num_slabs);
    subvol_shape = [vol_shape(1:2), num_sli];
    PrintVerbose(verbose, '\n rotation axis offset: %f', rot_axis_offset );
    PrintVerbose(verbose, '\n rotation axis position: %f', rot_axis_pos );    
    PrintVerbose(verbose, '\n max shape of reconstruction volume: [%g, %g, %g]', vol_shape )
    PrintVerbose(verbose, '\n shape of reconstruction subvolume: [%g, %g, %g]', subvol_shape )
    PrintVerbose(verbose, '\n available gpu memory : %g MiB (%.2f%%) of %g MiB', gpu.AvailableMemory/10^6, 100*gpu.AvailableMemory/gpu.TotalMemory, gpu.TotalMemory/10^6)
    PrintVerbose(verbose, '\n memory required for reconstruction of a single slice: %g MiB (estimate)', rec_mem / 1024^2 )
    PrintVerbose(verbose, '\n memory of reconstructed volume: %g MiB', prod( vol_shape ) * 4 / 1024^2 )
    PrintVerbose(verbose, '\n max memory of subvolume: %g MiB', prod( subvol_shape ) * 4 / 1024^2 )
    PrintVerbose(verbose, '\n number of subvolume slabs: %g', num_slabs )
    PrintVerbose(verbose, '\n max number of slices per slab: %g', num_sli )
        
    PrintVerbose(verbose, '\nTomographic reconstruction of %u slabs:', num_slabs)
    if isempty( take_neg_log )
        take_neg_log = 1;
    end
    vol_min = NaN;
    vol_max = NaN;
    angles_reco = angles;
    num_proj_reco = numel( angles );
    if isequal( angles(1), angles(end) )
        angles_reco(end) = [];
        num_proj_reco = numel( angles_reco );
    end
    for nn = 1:num_slabs

        % Slice indices
        sli0 = 1 + ( nn - 1) * num_sli;
        sli1 = min( [nn * num_sli, vol_shape(3)] );

        % Filter sinogram
        PrintVerbose(verbose, '\n slab %g of %g: Filter sino: ', nn, num_slabs)        
        filt = iradonDesignFilter(fbp_filter_type, (1 + fbp_filter_padding) * vol_shape1, fpb_freq_cutoff);       
        if butterworth_filter(1)
            [b, a] = butter(butterworth_order, butterworth_cutoff_frequ);
            bw = freqz(b, a, numel(filt) );
            filt = filt .* bw;
        end
        if fbp_filter_padding(1)
            sino = padarray( NegLog(permute(proj(:,sli0:sli1,1:num_proj_reco), [1 3 2]), take_neg_log), [vol_shape1 0 0], 'symmetric', 'post' );
            sino = real( ifft( bsxfun(@times, fft( sino, [], 1), filt), [], 1, 'symmetric') );
            sino = sino(1:vol_shape1,:,:);
        else
            sino = real( ifft( bsxfun(@times, fft( NegLog(permute(proj(:,sli0:sli1,1:num_proj_reco), [1 3 2]), take_neg_log), [], 1), filt), [], 1, 'symmetric') );
        end        
        PrintVerbose(verbose, 'done.')

        % Backprojection
        PrintVerbose(verbose, ' Backproject: ')    
        subvol_shape(3) = size( sino, 3);
        subvol_size = vol_size;
        subvol_size(5) = subvol_shape(3) / vol_shape(3) * vol_size(5);
        subvol_size(6) = subvol_shape(3) / vol_shape(3) * vol_size(6);
        vol = astra_parallel3D( sino, rot_angle_offset + angles_reco, rot_axis_offset_reco, subvol_shape, subvol_size, astra_pixel_size, link_data, rot_axis_tilt);
        PrintVerbose(verbose, 'done.')

        % Save subvolume
        if write_reco(1)
            PrintVerbose(verbose, ' Save slices:')
            if do_phase_retrieval(1)
                save_path = reco_phase_path;
            else
                save_path = reco_path;
            end
            CheckAndMakePath(save_path)
            
            % Save reco path and scan path in file
            filename = [userpath, filesep, 'experiments/p05/pathtolastreco'];
            fid = fopen( filename , 'w' );
            fprintf( fid, '%s', save_path );   
            fclose( fid );
            
            % Save slices
            parfor ii = 1:size( vol, 3)
                filename = sprintf( '%sreco_%06u.tif', save_path, sli0 + ii - 1);
                write32bitTIFfromSingle( filename, vol( :, :, ii) )
            end
            PrintVerbose(verbose, ' done.')
        end
        vol_min = min( min( vol(:) ), vol_min );
        vol_max = max( max( vol(:) ), vol_max );
                        
    end
    PrintVerbose(verbose, ' Time elapsed: %.1f s (%.2f min)', toc-t, (toc-t)/60 )    
    
    % Write log file
    if write_reco(1)
         logfile_name = sprintf( '%sreco.log', save_path);
         fid = fopen(logfile_name, 'w');
         fprintf(fid, 'scan_name : %s\n', scan_name);
         fprintf(fid, 'beamtime_id : %s\n', beamtime_id);
         fprintf(fid, 'scan_path : %s\n', scan_path);
         fprintf(fid, 'reco_path : %s\n', save_path);                  
         fprintf(fid, 'MATLAB notation, index of first element: 1, range: first:stride:last\n');
         fprintf(fid, 'MATLAB version : %s\n', version);
         fprintf(fid, 'platform : %s\n', computer);
         fprintf(fid, 'camera : %s\n', cam);
         fprintf(fid, 'num_dark_found : %u\n', num_dark);
         fprintf(fid, 'num_ref_found : %u\n', num_ref_found);
         fprintf(fid, 'num_ref_used : %u\n', num_ref_used);
         fprintf(fid, 'ref_range : %u:%u:%u\n', ref_range(1), ref_range(2) - ref_range(1), ref_range(end) );         
         fprintf(fid, 'num_proj_found : %u\n', num_proj_found);
         fprintf(fid, 'num_proj_used : %u\n', num_proj_used);
         fprintf(fid, 'proj_range : %u:%u:%u\n', proj_range(1), proj_range(2) - proj_range(1), proj_range(end) );
         fprintf(fid, 'raw_image_shape : %u %u\n', raw_im_shape);
         fprintf(fid, 'raw_image_shape_binned : %u %u\n', raw_im_shape_binned);
         fprintf(fid, 'binning_factor : %u\n', bin);
         fprintf(fid, 'effective_pixel_size_mu : %g\n', eff_pixel_size * 1e6);
         fprintf(fid, 'effective_pixel_size_binned_mu : %g\n', eff_pixel_size_binned * 1e6);
         fprintf(fid, 'energy : %g eV\n', energy);
         fprintf(fid, 'flat_field_correlation_area_1 : %u:%u:%u\n', flat_corr_area1(1), flat_corr_area1(2) - flat_corr_area1(1), flat_corr_area1(end));
         fprintf(fid, 'flat_field_correlation_area_2 : %u:%u:%u\n', flat_corr_area2(1), flat_corr_area2(2) - flat_corr_area2(1), flat_corr_area2(end));
         if ~read_flatcor
             fprintf(fid, 'min_max_of_all_darks : %6g %6g\n', dark_min, dark_max);
             fprintf(fid, 'min_max_of_median_dark : %6g %6g\n', dark_med_min, dark_med_max);
             fprintf(fid, 'min_max_of_all_flats : %6g %6g\n', flat_min, flat_max);
             fprintf(fid, 'min_max_of_all_corrected_flats : %6g %6g\n', flat_min2, flat_max2);
             fprintf(fid, 'min_max_of_all_raws :  %6g %6g\n', raw_min, raw_max);         
             fprintf(fid, 'min_max_of_all_corrected_raws :  %6g %6g\n', raw_min2, raw_max2);
             fprintf(fid, 'min_max_of_all_flat_corr_projs : %g %g \n', proj_min, proj_max);
         end
         % Phase retrieval
         fprintf(fid, 'do_phase_retrieval : %u\n', do_phase_retrieval);
         fprintf(fid, 'phase_retrieval_method : %s\n', phase_retrieval_method);
         fprintf(fid, 'phase_retrieval_regularisation_parameter : %f\n', phase_retrieval_reg_par);
         fprintf(fid, 'phase_retrieval_binary_filter_threshold : %f\n', phase_retrieval_bin_filt);
         fprintf(fid, 'phase_padding : %u\n', phase_padding);
         % Rotation
         fprintf(fid, 'rotation_angle_full_rad : %f\n', rot_angle_full);
         fprintf(fid, 'rotation_angle_offset_rad : %f\n', rot_angle_offset);
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
         fprintf(fid, 'full_reconstruction_time_s : %.1f\n', toc);
         fprintf(fid, 'date : %s', datetime);
         fclose(fid);
         PrintVerbose(verbose, '\n log file : %s', logfile_name)
    end
end

PrintVerbose(verbose, '\nFinished. Total time elapsed: %g s (%.2f min)\n\n', toc, toc / 60 );
% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
