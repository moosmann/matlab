function p05_reco(external_parameter)
% P05/P07 reconstruction pipeline: preprocessing,filtering,phase retrieval,
% tomographic reconstruction,...
%
% USAGE:
% Edit parameters in PARAMETERS / SETTINGS section below and run script.
%
% HOW TO RUN THE SCRIPT:
% - Editor windows: press 'F5' when focus is in the Editor window
% - Editor tab: click 'Run' in the toolstrip
% - Command Window: type 'p05_reco' and hit Enter
%
% HOW TO AUTOMATICALLY LOOP RECO OVER DATA SETS:a
% To loop over different data or parameters sets see
% 'p05_reco_loop_template' and/or 'p05_create_reco_loop_script'.
%
% For additional information see 'p05_reco_NOTES'.
%
% Please cite following article in the case of publication:
% - Moosmann,J. et al. Time-lapse X-ray phase-contrast microtomography for
% in vivo imaging and ansdfadsfaalysis of morphogenesis Nat. Protocols 9,
% 294-304 (2014)
% - Moosmann,J. moosmann/matlab:. Zenodo. https:// doi. org/ 10. 5281/ ZENODO. 51187 37 (2021).
% - ASTRA Toolbox,see http://www.astra-toolbox.com/
%
% Latest version: https://github.com/moosmann/matlab.git
%
% Written by Julian Moosmann.

dbstop if error

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETERS / SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% !!! QUICK SWITCH TO ALTERNATIVE SET OF PARAMETERS !!!
% !!! OVERWRITES PARAMETERS BELOW QUICK SWITCH SECTION !!!
% Just copy parameter and set quick switch to 1
par.quick_switch = 0;

% END OF QUICK SWITCH TO ALTERNATIVE SET OF PARAMETERS %%%%%%%%%%%%%%%%%%%%

pp_parameter_switch % DO NOT DELETE OR EDIT THIS LINE %%%%%%%%%%%%%%%%%%%%%

%%% SCAN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.scan_path = pwd; % string/pwd. pwd: change to directory of the scan to be reconstructed,string: absolute scan path,last_folder_modified('folder')
par.ref_path = {}; % cell of strings. Additonal data sets to be included for the correlation of projections and reference images
par.nexus_path = ''; % string,absolute path to h5 file
par.read_flatcor = 0; % read preprocessed flatfield-corrected projections. CHECK if negative log has to be taken!
par.read_flatcor_path = ''; % absolute path containing flat-field corrected projections
par.read_flatcor_range = 1; % scalar or vector. range of flatcorrected projections to be read
par.read_flatcor_bin = 1; % Binning of flat-corrected projections
par.read_flatcor_trafo = @(im) im; %fliplr(im); % anonymous function applied to the image which is read e.g. @(x) rot90(x)
par.read_sino = 0; % read preprocessed sinograms. CHECK if negative log has to be taken!
par.read_sino_folder = ''; % subfolder to scan path
par.read_sino_trafo = @(x) (x);%rot90(x); % anonymous function applied to the image which is read e.g. @(x) rot90(x)
par.read_sino_range = 1; % scalar,vector. If integer in [1,10]: interpreted as increment,else slice number. if scalar in [0,1): interpreted as relative slice number. if 2-vector: interpreted as absolute or relative slice number of first/last slice to read
par.sino_roi = []; % horizontal ROI when reading sinograms,vertical ROI not yet implemented
par.filter_sino = 1; % bool. Pixel filtering of sinogram using parameters below:
pixel_filter_sino.threshold_hot = 0;
pixel_filter_sino.threshold_dark = 0;
pixel_filter_sino.medfilt_neighboorhood = [3 3];
pixel_filter_sino.filter_dead_pixel = 1;
pixel_filter_sino.filter_Inf = 1;
pixel_filter_sino.filter_NaN = 1;
pixel_filter_sino.verbose = 0;
par.energy = []; % eV! if empty: read from log file (log file values can be ambiguous or even missing sometimes)
par.sample_detector_distance = []; % in m. if empty: read from log file
par.eff_pixel_size = []; % in m. if empty: read from log lfile. effective pixel size =  detector pixel size / magnification
par.pixel_scaling = []; % to account for mismatch of eff_pixel_size with,ONLY APPLIED BEFORE TOMOGRAPHIC RECONSTRUCTION,HAS TO BE CHANGED!
par.read_image_log = 0; % bool,default: 0. Read metadata from image log instead hdf5,if image log exists
par.read_filenames_from_disk = 0; % only for stepscans with tiff subfolders
%%% PREPROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.raw_bin = 2; % projection binning factor: integer
par.raw_roi = []; % vertical and/or horizontal ROI; coordinate (1,1) = top left pixel; supports absolute,relative,negative,and mixed indexing.
% []: use full image; oi=-1 defaults to min(proj(:,:,[1 end])) + 4*median(dark(:))
par.im_trafo = '';% 'rot90(im,1)'; % string to be evaluated after reading data in the case the image is flipped/rotated/etc due to changes at the beamline,e.g. 'rot90(im)'
par.filter_ref = @(x) (x);
par.filter_proj = @(x) (x);
% STITCHING/CROPPING only for scans without lateral movment. Legacy support
par.crop_at_rot_axis = 0; % for recos of scans with excentric rotation axis but WITHOUT projection stitching
par.stitch_projections = 0; % for 2 pi cans: stitch projection at rotation axis position. Recommended with phase retrieval to reduce artefacts. Standard absorption contrast data should work well without stitching. Subpixel stitching not supported (non-integer rotation axis position is rounded,less/no binning before reconstruction can be used to improve precision).
par.stitch_align_overlap = 0; %25; % pos. integer,number of pixels,if 0 do nothing,align mean values of vertical aread left/right to the rot axis within given pixel range. Use with care: can reduce, but also enhance artefacts significantly. Typically only works well for off-centered scans.
par.stitch_method = 'step';'linear';'sine'; %  ! CHECK correlation area !
% 'step' : no interpolation,use step function
% 'linear' : linear interpolation of overlap region
% 'sine' : sinusoidal interpolation of overlap region
par.proj_range = []; % range of projections to be used (from all found,if empty or 1: all,if scalar: stride,if range: start:incr:end
par.ref_range = [];% range of flat fields to be used (from all found),if empty or 1: all. if scalar: stride,if range: start:incr:end
par.crop_proj = 0; % Crop images to account for random lateral shift
par.virt_s_pos = 0; % Correct sample position in reconsructed volume if virtual sample position motors are used
pixel_filter_threshold_dark = [0.02 0.005]; % Dark fields: threshold parameter for hot/dark pixel filter,for details see 'FilterPixel'
pixel_filter_threshold_flat = [0.02 0.005]; % Flat fields: threshold parameter for hot/dark pixel filter,for details see 'FilterPixel'
pixel_filter_threshold_proj = [0.02 0.02]; % Raw projection: threshold parameter for hot/dark pixel filter,for details see 'FilterPixel'
pixel_filter_radius = [5 5]; % Increase only if blobs of zeros or other artefacts are expected. Can increase processing time heavily.
par.ring_current_normalization = 0; % normalize flat fields and projections by ring current
image_correlation.method = 'median';'ssim-ml';'entropy';'none';'ssim';'ssim-g';'std';'cov';'corr';'diff1-l1';'diff1-l2';'diff2-l1';'diff2-l2';'cross-entropy-12';'cross-entropy-21';'cross-entropy-x';
% Correlation of projections and flat fields. Essential for DCM data. Typically improves reconstruction quality of DMM data,too.
% Available methods ('ssim-ml'/'entropy' usually work best):
% 'none' : no correlation
% 'mean'/'median': mean/median flat
% 'ssim-ml' : Matlab's structural similarity index (SSIM),includes Gaussian smoothing
% 'entropy' : entropy measure of proj over flat,usually similar result to SSIM,but faster
% 'ssim' : own implementation of SSIM without smoothing,usually worse,but sometimes better than 'ssim-ml'
% 'ssim-g' : 'ssim' with smoothing (Gaussian blurring)
% 'cov' : cross covariance
% 'corr' : cross correlation = normalized cross covariance
% 'std' : standard deviation of proj over flat
% 'diff1/2-l1/2': L1/L2-norm of anisotropic (diff1-l*) or isotropic (diff2-l*) difference of projections and flat fields
% 'cross-entropy-*' : asymmetric (12,21) and symmetric (x) cross entropy
image_correlation.force_calc = 0; % bool. force compuation of correlation even though a (previously computed) corrlation matrix exists
image_correlation.num_flats = 11; % integer. number of best maching flat fields used for correction
image_correlation.area_width = [1 100];% 2-vector. correlation area: index vector or relative/absolute position of [first pix,last pix],negative indexing is supported
image_correlation.area_height = [0.25 0.75]; % 2/vector. correlation area [bottom top]: index vector or relative/absolute position of [first pix,last pix]
image_correlation.filter = 1; % bool,filter ROI before correlation
image_correlation.filter_type = 'median'; % string. correlation ROI filter type,currently only 'median' is implemnted
image_correlation.filter_parameter = {[5 5],'symmetric'}; % cell. filter paramaters to be parsed with {:}
% 'median' : using medfilt2,parameters: {[M N]-neighboorhood,'padding'}
% 'wiener' : using wiender2,parameters: {[M N]-neighboorhood}
ring_filter.apply = 1; % ring artifact filter (use only for scans without lateral sample movement)
ring_filter.apply_before_stitching = 0; % ! Consider when phase retrieval is applied !
ring_filter.method = 'jm'; %'wavelet-fft';'jm';
ring_filter.waveletfft_dec_levels = 1:6; % decomposition levels for 'wavelet-fft'
ring_filter.waveletfft_wname = 'db7';'db25';'db30'; % wavelet type,see 'FilterStripesCombinedWaveletFFT' or 'waveinfo'
ring_filter.waveletfft_sigma = 3; % integer scalar. suppression factor for 'wavelet-fft'
ring_filter.jm_median_width = 11; % integer scalar or vector. median averaging filter to be applied to angular averaged sinogram,multiple widths are applied consecutively,eg [3 11 21 31 39];
par.strong_abs_thresh = 1; % if 1: does nothing,if < 1: flat-corrected values below threshold are set to one. Try with algebratic reco techniques.
par.delete_empty_projections = 0; % bool. Delete projections that conain zeros
par.norm_sino = 0; % not recommended,can introduce severe artifacts,but sometimes improves quality
% Workaround correction for image distortions using a quadratic dilation/compression of the projections/sinogram
% Preferably,projection cropping of laterally shifted projection 'crop_proj' is not used
par.distortion_correction_distance = 0; % scalar,in binned pixel,distance between two regions in the tomogram that can be properly reconstructed using different rotation axis offsets,if 0: no correction done
par.distortion_correction_outer_offset = 0; % scalar,in pixel,rotation axis offset for the outer region. the offset for the inner region is used for reconstruction
par.distortion_correction_exponent = 2; % scalar, exponent of interpolation function: xq = x - 2 * offset_diff * (x / dist_offset).^exponent;
%%% PHASE RETRIEVAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phase_retrieval.apply = 0; % bool. See 'PhaseFilter' for detailed description of parameters !
phase_retrieval.apply_before = 0; % bool. before stitching,interactive mode,etc. For phase-contrast data with an excentric rotation axis phase retrieval should be done afterwards. To find the rotataion axis position use this option in a first run,and then turn it of afterwards.
phase_retrieval.post_binning_factor = 1; % bool. Binning factor after phase retrieval,but before tomographic reconstruction
phase_retrieval.method = 'tie'; % string. Available methods: 'qp' 'ctf' 'tie' 'qp2' 'qpcut','tieNLO_Schwinger'
% Interactive phase retrieval not supported for method 'tieNLO_Schwinger'
phase_retrieval.reg_par = 1.0; % regularization parameter. larger values tend to blurrier images. smaller values tend to original data.
phase_retrieval.bin_filt = 0.1; % threshold for quasiparticle retrieval 'qp','qp2'
phase_retrieval.cutoff_frequ = 2 * pi; % in radian. frequency cutoff in Fourier space for 'qpcut' phase retrieval
phase_retrieval.padding = 1; % padding of intensities before phase retrieval,0: no padding
phase_retrieval.tieNLO_Schwinger.sn = 10; % Schwinger regularization: points of support
phase_retrieval.tieNLO_Schwinger.smax = 10; % Schwinger regularization: maximumg support range
%%% TOMOGRAPHY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tomo.run = 1; % run tomographic reconstruction
tomo.run_interactive_mode = 1; % if tomo.run = 0,use to determine rot axis positions without processing the full tomogram;
tomo.reco_mode = '3D';'slice'; % string. slice-wise or full 3D backprojection. 'slice': volume must be centered at origin & no support of rotation axis tilt,reco binning,save compressed
tomo.vol_size = [];%[-1 1 -1 1 -0.5 0.5];% 6-component vector [xmin xmax ymin ymax zmin zmax] for excentric rot axis pos or extended FoV;. if empty,volume is centerd within tomo.vol_shape. unit voxel size is assumed. if smaller than 10 values are interpreted as relative size w.r.t. the detector size. Take care bout minus signs! Note that if empty vol_size is dependent on the rotation axis position.
% Orientation using Matlab's matrix notation: relative coordinates in a horizontal reconstruction plane using imagej (-0.5,-0.5) = top left
% pixel,(0.5,0.5) = bottom right,(0,0) = center
% [left,right,top,bottom,bottom slice,top slice]
tomo.vol_shape = []; %[1 1 1] integer vector. shape (# voxels) of reconstruction volume. used for excentric rot axis pos. if empty,inferred from 'tomo.vol_size'. in absolute numbers of voxels or in relative number w.r.t. the default volume which is given by the detector width and height.
tomo.rot_angle_full_range = [];% 2 * pi ;[]; % in radians. if []: full angle of rotation including additional increment,or array of angles. if empty full rotation angles is determined automatically to pi or 2 pi
tomo.rot_angle_offset = 0; % global rotation of reconstructed volume
tomo.interpolate_missing_angles = 0; % limited or missing angle tomography
tomo.rot_axis_offset = [] / 1 * par.raw_bin; % rotation axis offset w.r.t to the image center. Assuming the rotation axis position to be centered in the FOV for standard scan,the offset should be close to zero.
tomo.rot_axis_position = []; % if empty use automatic computation. EITHER OFFSET OR POSITION MUST BE EMPTY. YOU MUST NOT USE BOTH!
tomo.rot_axis_offset_shift = []; % absolute lateral movement in pixels during fly-shift-scan,overwrite lateral shift read out from hdf5 log
tomo.vert_shift = []; % vertical shift for spiral/helical CT
tomo.flip_scan_position = 0; % for debugging
tomo.rot_axis_tilt_camera = 0; % in rad. camera tilt w.r.t rotation axis.
tomo.rot_axis_tilt_lamino = 0; % in rad. lamino tilt w.r.t beam.
tomo.rot_axis_corr_area1 = []; % ROI to correlate projections at angles 0 & pi. Use [0.75 1] or so for scans with an excentric rotation axis
tomo.rot_axis_corr_area2 = []; % ROI to correlate projections at angles 0 & pi
tomo.fbp_filter_type = 'Ram-Lak';'linear'; % see iradonDesignFilter for more options. Ram-Lak according to Kak/Slaney
tomo.fbp_filter_freq_cutoff = 1; % Cut-off frequency in Fourier space of the above FBP filter
tomo.fbp_filter_padding = 1; % symmetric padding for consistent boundary conditions,0: no padding
tomo.fbp_filter_padding_method = 'symmetric';
tomo.butterworth_filter = 0; % use butterworth filter in addition to FBP filter
tomo.butterworth_filter_order = 1;
tomo.butterworth_filter_frequ_cutoff = 0.9;
tomo.astra_pixel_size = 1; % detector pixel size for reconstruction: if different from one 'tomo.vol_size' must to be ajusted,too!
tomo.take_neg_log = []; % take negative logarithm. if empty,use 1 for attenuation contrast,0 for phase contrast
tomo.algorithm =  'fbp';'cgls';'sirt';'sart';'em';'fbp-astra'; % SART/EM only work for 3D reco mode
tomo.iterations = 50; % for iterateive algorithms: 'sirt','cgls','sart','em'
tomo.MinConstraint = []; % sirt3D/sirt2d/sart2d only. If specified,all values below MinConstraint will be set to MinConstraint. This can be used to enforce non-negative reconstructions,for example.
tomo.MaxConstraint = []; % sirt3D/sirt2d/sart2d only. If specified,all values above MaxConstraint will be set to MaxConstraint.
tomo.rot_axis_search_auto = 0; % find extrema of metric within search range
tomo.rot_axis_search_range = []; % search reach for automatic determination of the rotation axis offset,overwrite interactive result if not empty
tomo.rot_axis_search_metric = 'iso-grad'; % string: 'neg','entropy','iso-grad','laplacian','entropy-ML','abs'. Metric to find rotation axis offset
tomo.rot_axis_search_extrema = 'max'; % string: 'min'/'max'. chose min or maximum position
tomo.rot_axis_search_fit = 1; % bool: fit calculated metrics and find extrema,otherwise use extrema from search range
tomo.rot_axis_offset_metric_roi = []; % 4-vector: [. ROI for metric calculation. roi = [y0,x0,y1-y0,x1-x0]. (x,y)=(0,0)=upper left
tomo.rot_axis_search_slice = []; % scalar: slice used to find rot axis. if empty: uses slice from interactive mode,if that is empty uses central slice.
tomo.rot_axis_search_range_from_interactive = 0; % boolean: use search range from interactive mode
%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write.path = ''; %'/gpfs/petra3/scratch/moosmanj'; % string. absolute path were output data will be stored. !!overwrites the write.to_scratch flag. if empty uses the beamtime directory and either 'processed' or 'scratch_cc'
write.to_scratch = 0; % write to 'scratch_cc' instead of 'processed'
write.deleteFiles = 0; % delete files already existing in output folders. Useful if number or names of files differ when reprocessing.
write.beamtimeID = ''; % string (regexp),typically beamtime ID,mandatory if 'write.deleteFiles' is true (safety check)
write.scan_name_appendix = ''; % appendix to the output folder name which defaults to the scan name
write.parfolder = ''; % parent folder to 'reco','sino','phase',and 'flat_corrected'
write.subfolder_flatcor = ''; % subfolder in 'flat_corrected'
write.subfolder_phase_map = ''; % subfolder in 'phase_map'
write.subfolder_sino = ''; % subfolder in 'sino'
write.subfolder_reco = ''; % subfolder in 'reco'
write.flatcor = 0; % save preprocessed flat corrected projections
write.flatcor_stitched = 0; % save stitched flat corrected projections (when stitchting)
write.phase_map = 0; % save phase maps (if phase retrieval is not 0)
write.sino = 0; % save sinograms (after preprocessing & before FBP filtering and phase retrieval)
write.phase_sino = 0; % save sinograms of phase maps
write.reco = 1; % save reconstructed slices (if tomo.run=1)
write.float = 1; % single precision (32-bit float) tiff
write.float_adapthisteq = 0; % save float with adaptive histogram equalization filter. only for 3D recos currently
write.uint16 = 0; % save 16bit unsigned integer tiff using 'write.compression_method'
write.uint8 = 0; % save binned 8bit unsigned integer tiff using 'write.compression_method'
% Optionally save binned reconstructions,only works in '3D' reco_mode
write.float_binned = 0; % save binned single precision (32-bit float) tiff
write.uint16_binned = 0; % save binned 16bit unsigned integer tiff using 'write.comression_method'
write.uint8_binned = 0; % save binned 8bit unsigned integer tiff using 'wwrite.compression_method'
write.reco_binning_factor = 2; % IF BINNED VOLUMES ARE SAVED: binning factor of reconstructed volume
write.compression_method = 'outlier';'histo';'full'; 'std'; 'threshold'; % method to compress dynamic range into [0,1]
write.compression_parameter = [0.02 0.02]; % compression-method specific parameters
% methods for the compression of the dynamic range:
% 'outlier' : [LOW,HIGH] = write.compression_parameter,eg. [0.01 0.03],outlier given in percent,if scalar LOW = HIGH.
% 'full' : full dynamic range is used
% 'threshold' : [LOW HIGH] = write.compression_parameter,eg. [-0.01 1]
% 'std' : NUM = write.compression_parameter,mean +/- NUM*std,dynamic range is rescaled to within -/+ NUM standard deviations around the mean value
% 'histo' : [LOW HIGH] = write.compression_parameter (100*LOW)% and (100*HIGH)% of the original histogram,e.g. [0.02 0.02]
write.uint8_segmented = 0; % experimental: threshold segmentaion for histograms with 2 distinct peaks: __/\_/\__
write.outputformat = 'tif';'hdf_volume'; % string. Not yet implemented for all reco modes
%%% INTERACTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.visual_output = 1; % show images and plots during reconstruction
interactive_mode.rot_axis_pos = 1; % reconstruct slices with dif+ferent rotation axis offsets
interactive_mode.rot_axis_pos_default_search_range = []; % if empty: asks for search range when entering interactive mode
interactive_mode.rot_axis_tilt = 0; % reconstruct slices with different offset AND tilts of the rotation axis
interactive_mode.rot_axis_tilt_default_search_range = []; % if empty: asks for search range when entering interactive mode
interactive_mode.lamino = 1; % find laminography tilt instead camera tilt
interactive_mode.angles = 0; % reconstruct slices with different scalings of angles
interactive_mode.angle_scaling_default_search_range = []; % if empty: use a variaton of -/+5 * (angle increment / maximum angle)
interactive_mode.slice_number = 0.5; % default slice number. if in [0,1): relative,if in (1,N]: absolute
interactive_mode.phase_retrieval = 1; % Interactive retrieval to determine regularization parameter
interactive_mode.phase_retrieval_default_search_range = []; % if empty: asks for search range when entering interactive mode,otherwise directly start with given search range
interactive_mode.show_stack_imagej = 1; % use imagej instead of MATLAB to scroll through images during interactive mode
interactive_mode.show_stack_imagej_use_virtual = 0; % use virtual stack for faster loading,but slower scrolling
%%% HARDWARE / SOFTWARE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tomo.astra_link_data = 0; % boolean: ASTRA data objects become references to Matlab arrays. Reduces memory issues.
par.gpu_index = []; % integer vector: indices of GPU devices to use,Matlab notation: index starts from 1. default: [],uses all
par.use_cluster = 0; % if available: on MAXWELL nodes disp/nova/wga/wgs cluster computation can be used. Recommended only for large data sets since parpool creation and data transfer implies a lot of overhead.
par.use_gpu_in_parfor = 0; % boolean
pixel_filter_sino.use_gpu = par.use_gpu_in_parfor;
par.poolsize = 50; % scalar: number of workers used in a local parallel pool. if 0: use current config. if >= 1: absolute number. if 0 < poolsize < 1: relative amount of all cores to be used. if SLURM scheduling is available,a default number of workers is used.
par.poolsize_gpu_limit_factor = 0.5; % scalar: relative amount of GPU memory used for preprocessing during parloop. High values speed up Proprocessing,but increases out-of-memory failure
phase_retrieval.use_parpool = 1; % bool. Disable parpool when out-of-memory error occurs during phase retrieval.
par.window_state = 'minimized';'normal';'maximized';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% END OF PARAMETERS / SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('START RECONSTRUCTION: ')

weblink_url = 'https://github.com/moosmann/matlab';
weblink_name = sprintf('%s',weblink_url);
weblink = sprintf('<a href = "%s">%s</a>\n',weblink_url,weblink_name);
fprintf('\n code repository on github: ')
fprintf(weblink)

close all hidden;
close all force;

%%% Parameters set by reconstruction loop script 'p05_reco_loop' %%%%%%%%%%
if nargin == 1 %exist('external_parameter','var')
    % Fields of parameter struct from loop script
    field_name_cell = fieldnames(external_parameter);
    for nn = 1:numel(field_name_cell)
        field_name = field_name_cell{nn};
        field_value = external_parameter.(field_name);
        eval(sprintf('%s = field_value;',field_name));
    end
    %clear external_parameter field_name_cell field_name field_value
    par.quick_switch = 0;
    %par.visual_output = 1;
end

%%% QUICK SWITCH PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(par,'quick_switch') && par.quick_switch
    % Loop over parameter structs
    for sn = 1:par_quick_switch.num_structs
        % Parameter struct name
        struct_name = par_quick_switch.who{sn};
        if strcmp(struct_name,'ans')
            continue
        end
        % Parameter struct
        struct_tmp = par_quick_switch.structs.(struct_name);
        % Fields of parameter struct
        field_name_cell = fieldnames(struct_tmp);
        % Loop over parameter struct fields
        for fn = 1:numel(field_name_cell)
            % current parameter struct field
            field_name = field_name_cell{fn};
            % Evaluate parameter.fields
            eval_string = sprintf('%s.(''%s'') = par_quick_switch.structs.(''%s'').(''%s'');',struct_name,field_name,struct_name,field_name);
            eval(eval_string);
        end
    end
    clear par_quick_switch struct_name struct_tmp field_name_cell field_name eval_string
    cprintf('Red','\nATTENTION: Quick parameter switch is turned on!\n\n')
end
tic;
verbose = 1;
vert_shift = 0;
offset_shift = 0;
scan_position = [];
logpar = [];
s_stage_z_str =  '';
s_stage_z_setup = [];
s_stage_x_setup = [];
imlogcell = [];
scan_position_index = [];
par.raw_data = 0;
par.verbose = verbose;
par.s_in_pos_mm = 0;

% Parameter checks
if ~phase_retrieval.apply
    phase_retrieval.post_binning_factor = 0;
end
if interactive_mode.rot_axis_tilt && strcmpi(tomo.reco_mode,'slice')
    error('Slicewise reconstruction and reconstruction with tilted rotation axis are not compatible!')
end
% infer (relative) vol_shape from vol_size if vol_shape is empty,but vol_size is given
if ~isempty(tomo.vol_size) && isempty(tomo.vol_shape)
    tomo.vol_shape = tomo.vol_size(2:2:end) - tomo.vol_size(1:2:end);
end
if ~isempty(tomo.rot_axis_offset) && ~isempty(tomo.rot_axis_position)
    error('tomo.rot_axis_offset (%f) and tomo.rot_axis_position (%f) cannot be used simultaneously. One must be empty.',tomo.rot_axis_offset,tomo.rot_axis_position)
end
% Take negative logarithm for attenuation contrast,unless phase contrast
% is used
if isempty(tomo.take_neg_log)
    if phase_retrieval.apply
        tomo.take_neg_log = 0;
    else
        tomo.take_neg_log = 1;
    end
end

%% Default assignment if non-existing or empty!
%assign_default('',)
assign_default('par.raw_bin',2);
assign_default('par.read_flatcor_bin',1);
assign_default('par.read_flatcor_range',1);
assign_default('par.im_shape_raw',[]);
assign_default('par.im_format','');
assign_default('par.tif_info',[]);
assign_default('par.im_trafo','');
assign_default('par.dtype','');
assign_default('par.energy',[]);
assign_default('par.sample_detector_distance',[]);
assign_default('par.eff_pixel_size',[]);
assign_default('par.pixel_scaling',[]);
assign_default('par.ref_range',1)
assign_default('par.proj_range',1)
assign_default('par.read_image_log',0)
assign_default('par.read_filenames_from_disk',0)
assign_default('par.crop_at_rot_axis',0);
assign_default('par.stitch_projections',0);
assign_default('par.stitch_method','step');
assign_default('par.stitch_align_overlap',25);
assign_default('par.virt_s_pos',0);
assign_default('par.crop_proj',0);
assign_default('par.ring_current_normalization',1);
assign_default('pixel_filter_radius',[3 3])
assign_default('pixel_filter_sino',[]);
assign_default('image_correlation.filter',1);
assign_default('image_correlation.filter_type','median');
assign_default('image_correlation.filter_parameter',{[3 3],'symmetric'});
assign_default('image_correlation.force_calc',0);
assign_default('phase_retrieval.use_parpool',0);
assign_default('tomo.reco_mode','3D')
assign_default('tomo.rot_axis_offset',0)
assign_default('tomo.rot_axis_tilt_camera',0)
assign_default('tomo.rot_axis_tilt_lamino',0)
assign_default('tomo.flip_scan_position',0)
assign_default('tomo.rot_axis_corr_area2',[0.1 0.9]);
assign_default('tomo.angle_scaling',1);
assign_default('tomo.MinConstraint',[])
assign_default('tomo.MaxConstraint',[])
assign_default('tomo.interpolate_missing_angles',0)
assign_default('tomo.vert_shift',[])
assign_default('tomo.rot_axis_search_auto',0);
assign_default('tomo.rot_axis_search_range',[]);
assign_default('tomo.rot_axis_search_metric','');
assign_default('tomo.rot_axis_search_verbose',1);
assign_default('tomo.rot_axis_search_extrema','max');
assign_default('tomo.rot_axis_search_fit',1);
assign_default('tomo.rot_axis_offset_metric_roi',[]);
assign_default('tomo.rot_axis_search_slice',[]);
assign_default('tomo.rot_axis_search_range_from_interactive',0);
assign_default('tomo.interactive_offset_range',[]);
assign_default('tomo.slice',[]);
assign_default('write.path','')
assign_default('write.flatcor',0)
assign_default('write.flatcor_stitched',0)
assign_default('write.parfolder','')
assign_default('write.subfolder_reco','')
assign_default('write.subfolder_flatcor','')
assign_default('write.subfolder_phase_map','')
assign_default('write.subfolder_sino','')
assign_default('write.sino_shift_cropped',0)
assign_default('write.deleteFiles',0)
assign_default('write.beamtimeID','')
assign_default('write.scan_name_appendix','')
assign_default('write.uint8_segmented',0)
assign_default('write.phase_appendix','')
assign_default('write.float_adapthisteq',1)
assign_default('write.outputformat','tif')
assign_default('interactive_mode.rot_axis_pos_default_search_range',-4:0.5:4)
assign_default('interactive_mode.rot_axis_tilt_default_search_range',-0.005:0.001:0.005)
assign_default('interactive_mode.rot_axis_search_range',[])
assign_default('interactive_mode.phase_retrieval',0)
assign_default('interactive_mode.phase_retrieval_default_search_range',[])
assign_default('interactive_mode.angles',0);
assign_default('interactive_mode.angle_scaling_default_search_range',[]);
assign_default('interactive_mode.show_stack_imagej',1)
assign_default('interactive_mode.show_stack_imagej_use_virtual',1)
assign_default('par.distortion_correction_distance',0)
assign_default('par.distortion_correction_outer_offset',0)
assign_default('par.distortion_correction_exponent',2)
assign_default('par.gpu_index',[])
assign_default('par.sino_roi',[])
assign_default('par.read_sino_range',1);
assign_default('par.nexus_path','');
assign_default('par.delete_empty_projections',0)
assign_default('par.window_state','normal');
assign_default('par.strong_abs_thresh',1);
assign_default('par.norm_sino',0);
%assign_default('',)

% Define variables from struct fields for convenience
par.raw_bin = single(par.raw_bin);
if ~par.read_flatcor
    raw_bin = par.raw_bin;
else
    raw_bin = par.read_flatcor_bin;
    par.raw_bin = raw_bin;
    if ~isempty(par.eff_pixel_size)
        if ~isempty(par.pixel_scaling)
            par.eff_pixel_size = par.pixel_scaling * par.eff_pixel_size;
        end
        par.eff_pixel_size_binned = raw_bin * par.eff_pixel_size;
    end
end
phase_bin = phase_retrieval.post_binning_factor;
window_state = par.window_state;
outputformat = write.outputformat;
exposure_time = [];
if isempty(par.energy)
    energy_was_empty = 1;
else
    energy_was_empty = 0;
end
cur = [];
astra_clear % if reco was aborted,ASTRA memory was not cleared
% Utility functions
imsc1 = @(im) imsc(rot90(im));
NameCellToMat = @(name_cell) reshape(cell2mat(name_cell),[numel(name_cell{1}),numel(name_cell)])';
% Disable warnings
warning('off','MATLAB:imagesci:rtifc:missingPhotometricTag');
warning('off','MATLAB:hg:AutoSoftwareOpenGL');
warning('off','parallel:gpu:device:DeviceDeprecated')
%% Folders
if par.read_flatcor
    %     while ~strcmp(b,'processed')
    %         [a,b]=fileparts(a);
    %     end
    par.scan_path = par.read_flatcor_path;
    write.flatcor = 0;
    cprintf('Red','\nReading flat corrected projections!\n')
end
% Scan path
if ~iscell(par.scan_path)
    while par.scan_path(end) == filesep
        par.scan_path(end) = [];
    end
    [par.raw_path,scan_name] = fileparts(par.scan_path);
    par.scan_path = [par.scan_path,filesep];
    scan_path = par.scan_path;
else
    a = par.scan_path{1};
    b = par.scan_path{2};
    while ~strcmp(a,b)
        a(end) = [];
        b(end) = [];
    end
    while a(end) == '_'
        a(end) = [];
    end
    [par.raw_path,scan_name] = fileparts(a);
    scan_path = [a,filesep];
    for n = 1:numel(par.scan_path)
        if par.scan_path{n}(end) ~= filesep
            par.scan_path{n} = [par.scan_path{n},filesep];
        end
    end
end
[beamtime_path,raw_folder] = fileparts(par.raw_path);
[~,beamtime_id] = fileparts(beamtime_path);
if ~strcmp(raw_folder,'raw') && ~par.read_sino && ~par.read_flatcor
    error('Given path does not contain a ''raw'' folder: %s',raw_folder)
end

% Output path and folder
out_folder = 'processed';
if write.to_scratch
    out_folder = 'scratch_cc';
end
if isempty(write.path)
    write.path = [beamtime_path,filesep,out_folder,filesep,scan_name];
else
    write.path = [write.path,filesep,scan_name];
end
if ~isempty(write.scan_name_appendix)
    write.path = [ write.path write.scan_name_appendix ];
end
write.parpath = [write.path filesep ];
if ~isempty(write.parfolder)
    write.path = [write.path,filesep,write.parfolder];
end
CheckAndMakePath(write.path);
fn_diary = sprintf('%s/command_window_diary.txt',write.path);
diary(fn_diary)
diary on

% Save raw path to file for shell short cut
filename = [userpath,filesep,'path_to_raw'];
%filename = [getenv('HOME'),filesep,'path_to_raw'];
fid = fopen(filename,'w');
fprintf(fid,'%s',par.raw_path);
fclose(fid);

fprintf('%s',scan_name)
write.scan_name = scan_name;
write.is_phase = phase_retrieval.apply;
fprintf(' at %s',datetime)
fprintf('\n scan_path:\n  %s',scan_path)
fprintf('\n provided nexus path: %s',par.nexus_path)

% Save scan path to file
filename = [userpath,filesep,'path_to_scan'];
%filename = [getenv('HOME'),filesep,'path_to_scan'];
fid = fopen(filename,'w');
fprintf(fid,'%s',scan_path);
fclose(fid);

% Reco path
if isempty(write.subfolder_reco)
    write.reco_path = [write.path,filesep,'reco',filesep];
else
    write.reco_path = [write.path,filesep,'reco',filesep,write.subfolder_reco,filesep];
end

% Figure path
write.fig_path = [write.path,filesep,'figures',filesep];
fig_path = write.fig_path;

% Image path
write.im_path = [write.path,filesep,'images',filesep];
im_path = write.im_path;
im_path1 = [write.im_path 'images1/'];
im_path2 = [write.im_path 'images2/'];
im_path3 = [write.im_path 'images3/'];
if phase_retrieval.apply
    im_path1reco = [write.im_path 'reco_phase1/'];
    im_path2reco = [write.im_path 'reco_phase2/'];
    im_path3reco = [write.im_path 'reco_phase3/'];
else
    im_path1reco = [write.im_path 'reco1/'];
    im_path2reco = [write.im_path 'reco2/'];
    im_path3reco = [write.im_path 'reco3/'];
end

% Path to flat-field corrected projections
flatcor_path = sprintf('%s/flat_corrected/rawBin%u/',write.path,raw_bin);
if par.stitch_projections == 1
    flatcor_path_stitched = sprintf('%s/flat_corrected/rawBin%u/',write.path,raw_bin);
end
if ~isempty(write.subfolder_flatcor)
    flatcor_path =  sprintf('%s%s/',flatcor_path,write.subfolder_flatcor);
    if par.stitch_projections == 1
        flatcor_path_stitched =  sprintf('%s%s/',flatcor_path_stitched,write.subfolder_flatcor);
    end
end
PrintVerbose(write.flatcor,'\n flatcor_path:\n  %s',flatcor_path)
% Path to retrieved phase maps
write.phase_map_path = sprintf('%s/phase_map/rawBin%u/',write.path,raw_bin);
if ~isempty(write.subfolder_phase_map)
    write.phase_map_path = sprintf('%s%s/', write.phase_map_path,write.subfolder_phase_map);
end
PrintVerbose(phase_retrieval.apply && write.phase_map,'\n phase_map_path:\n  %s',write.phase_map_path)
% Sinogram path
write.sino_par_path = sprintf('%s/sino/',write.path);
write.sino_path = sprintf('%s/sino/rawBin%u/',write.path,raw_bin);
write.sino_phase_path = sprintf('%s/sino_phase/rawBin%u/',write.path,raw_bin);
if ~isempty(write.subfolder_sino)
    write.sino_path = sprintf('%s%s/',write.sino_path,write.subfolder_sino);
    write.sino_phase_path = sprintf('%s%s/',write.sino_phase_path,write.subfolder_sino);
end
PrintVerbose(write.sino,'\n sino_path:\n  %s',write.sino_path)
PrintVerbose(phase_retrieval.apply & write.phase_sino,'\n sino_phase_path:\n  %s',write.sino_phase_path)

CheckAndMakePath(write.path)
filename = sprintf('%s/parameters.mat',write.path);
save(filename)
CheckAndMakePath(fig_path)
% interactive path
p = [write.path,filesep,'interactive',filesep];
p = regexprep(p,'processed','scratch_cc');
write.interactive_path = p;
% Memory
fprintf('\n user :  %s',getenv('USER'));
fprintf('\n hostname : %s',getenv('HOSTNAME'));
[mem_free,mem_avail_cpu,mem_total_cpu] = free_memory;
fprintf('\n RAM: free,available,total : %.0f GiB (%g%%),%.0f GiB (%g%%),%.0f GiB',round([mem_free/1024^3,100 * mem_free/mem_total_cpu,mem_avail_cpu/1024^3,100*mem_avail_cpu/mem_total_cpu mem_total_cpu/1024^3]))
% Start parallel CPU pool
par.pool_tmp_folder = [beamtime_path filesep 'scratch_cc'];
[~,par.poolsize] = OpenParpool(par.poolsize,par.use_cluster,par.pool_tmp_folder,0,[]);
if isempty(par.gpu_index)
    par.gpu_index = 1:gpuDeviceCount;
end
clearGPUs
% GPU info quick
fprintf('\n GPUs : [index,total memory/GiB] =\n ')
for mm = 1:numel(par.gpu_index)
    nn = par.gpu_index(mm);
    gpu = parallel.gpu.GPUDevice.getDevice(nn);
    mem_total_gpu = gpu.TotalMemory/1024^3;
    fprintf(' [%u %.3g]',nn,mem_total_gpu)
end
% GPU info detailed,needed for parpool optimization
mem_avail_gpu = zeros([1,numel(par.gpu_index)]);
mem_total_gpu = mem_avail_gpu;
ngpu = mem_avail_gpu;
parfor mm = 1:numel(par.gpu_index)
    nn = par.gpu_index(mm);
    ngpu(mm) = nn;
    gpu = gpuDevice(nn);
    gpu.reset;
    mem_avail_gpu(mm) = gpu.AvailableMemory;
    mem_total_gpu(mm) = gpu.TotalMemory;
end
for n = 1:numel(ngpu)
    ma = mem_avail_gpu(n)/1024^3;
    mt = mem_total_gpu(n)/1024^3;
    r = 100 * ma / mt;
    fprintf('\n GPU %u: memory: total: %.3g GiB,available: %.3g GiB (%.2f%%)',ngpu(n),mt,ma,r)
end
tomo.astra_gpu_index = par.gpu_index;
par.mem_avail_gpu = mem_avail_gpu;
par.mem_total_gpu = mem_total_gpu;

% Renderer
% r = rendererinfo;
% fprintf('\nGraphics renderer information')
% fprintf('\n GraphicsRenderer: %s',r.GraphicsRenderer)
% fprintf('\n RendererDevice: %s',r.RendererDevice)
% fprintf('\n Version: %s',r.Version)
% fprintf('\n HardwareSupportLevel: %s',r.Details.HardwareSupportLevel)
% if d.Software
%     fprintf('\n')
%     warning(' Software rendering is used. For improved GUI performance,log in directly to the Maxwell node with FastX to enable hardware acceleration.')
% end

if ~par.read_flatcor && ~par.read_sino
    % Projection range to read
    if isempty(par.proj_range)
        par.proj_range = 1;
    end
    %% Read image log
    imlog = dir(sprintf('%s*image.log',scan_path));
    if ~par.read_filenames_from_disk % default
        if par.read_image_log && ~isempty(imlog)
            imlog = [imlog.folder filesep imlog.name];
            fid = fopen(imlog);
            % name time image_key angle s_stage_x piezo petra
            imlogcell = textscan(fid,'%s%u64%u%f%f%f%f');%,'Delimiter',{'\n','\r'})
            fclose(fid);

            fns = imlogcell{1};
            im_time = imlogcell{2};
            im_key = imlogcell{3};
            im_angle = imlogcell{4} / 180 * pi;
            im_s_stage_x = imlogcell{5};
            %im_piezo = imlogcell{6};
            im_petra = imlogcell{7};

            stimg_name.scan.value = imlogcell{1};
            stimg_name.scan.time =  im_time;
            stimg_key.scan.value = imlogcell{3};
            stimg_key.scan.time = im_time;

            % PETRA ring current
            %petra.time = im_time(im_key==0);
            %petra.current = im_petra(im_key==0);
            petra.time = im_time;
            petra.current = im_petra;

            % rotation axis
            s_rot.time = im_time(im_key==0);
            s_rot.value = im_angle(im_key==0);

            % Read out lateral rotation axis shift form log file
            s_stage_x.time = im_time(im_key==0);
            s_stage_x.value = im_s_stage_x(im_key==0);

            offset_shift_mm = s_stage_x.value;
            % File names
            proj_names = fns(im_key==0)';
            ref_names = fns(im_key==1)';
            ref_full_path = cellfun(@(a) [scan_path a],ref_names,'UniformOutput',0);
            dark_names = fns(im_key == 2)';

            if isscalar(par.proj_range)
                par.num_ref_found = numel(proj_names);
                par.proj_range = 1:par.proj_range:par.num_proj_found;
            end
            % Angles %%
            angles = im_angle(im_key == 0);
            angles = angles(par.proj_range);
        else
            % note that if par.ref_path is not empty,it will may be changed
            %t = tic;
            %fprintf('\n Read file names from disk:')
            [proj_names,proj_full_path, ref_names,ref_full_path, dark_names,dark_full_path, par] =  pp_get_filenames(par);
            %fprintf(' %.1f s',toc - t)
        end
        % hdf5 log
        nexuslog_name = pp_get_nexuslog_names(par);
    else
        fn = dir([par.scan_path filesep 'tiff00*/*dar.tif']);
        for nn = numel(fn):-1:1
            [~,pf] = fileparts(fn(nn).folder);
            dark_names{nn} = [pf filesep fn(nn).name];
            dark_full_path{nn} = [fn(nn).folder filesep fn(nn).name];
        end
        fn = dir([par.scan_path filesep 'tiff00*/*ref.tif']);
        for nn = numel(fn):-1:1
            [~,pf] = fileparts(fn(nn).folder);
            ref_names{nn} = [pf filesep fn(nn).name];
            ref_full_path{nn} = [fn(nn).folder filesep fn(nn).name];
        end
        fn = dir([par.scan_path filesep 'tiff00*/*img.tif']);
        for nn = numel(fn):-1:1
            [~,pf] = fileparts(fn(nn).folder);
            proj_names{nn} = [pf filesep fn(nn).name];
            proj_full_path{nn} = [fn(nn).folder filesep fn(nn).name];
        end
        nexuslog_name = {''};
    end
    par.num_dark = numel(dark_names);
    par.num_ref_found = numel(ref_names);
    par.num_proj_found = numel(proj_names);
    % Ref range to read
    if isempty(par.ref_range)
        par.ref_range = 1;
    end
    if isscalar(par.ref_range)
        par.ref_range = 1:par.ref_range:par.num_ref_found;
    end
    if isscalar(par.proj_range)
        par.proj_range = 1:par.proj_range:par.num_proj_found;
    end
    %% TODO: fix for multi scan reco
    proj_full_path = proj_full_path(par.proj_range);
    ref_full_path = ref_full_path(par.ref_range);
    %dark_nums = CellString2Vec(dark_names);
    proj_nums = CellString2Vec(proj_names(par.proj_range));
    par.num_ref_used = numel(par.ref_range);
    par.num_proj_used = numel(par.proj_range);
    fprintf('\n refs found : %g',par.num_ref_found)
    fprintf('\n refs used : %g',par.num_ref_used)
    fprintf('\n reference range used : %g:%g:%g%',par.ref_range(1),par.ref_range(2) - par.ref_range(1),par.ref_range(end))
    fprintf('\n darks found : %g',par.num_dark)
    fprintf('\n projections found : %g',par.num_proj_found)
    fprintf('\n projections used : %g',par.num_proj_used)
    fprintf('\n projection range used : first:stride:last =  %g:%g:%g',par.proj_range(1),par.proj_range(2) - par.proj_range(1),par.proj_range(end))

    %% Scan Log %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % scanlog_name = sprintf('%s%sscan.log',scan_path,scan_name);
    % if ~exist(scanlog_name,'file')
    %     % back up for older log file name schemes
    %     str = dir(sprintf('%s*scan.log',scan_path));
    %     scanlog_name = sprintf('%s/%s',str.folder,str.name);
    % end
    % [logpar,cur,cam] = p05_log(scanlog_name);
    % if isempty(par.eff_pixel_size)
    %     par.eff_pixel_size = logpar.eff_pixel_size;
    % end
    % par.eff_pixel_size_binned = raw_bin * par.eff_pixel_size;
    % exposure_time = logpar.exposure_time;
    % if isfield(logpar,'s_in_pos')
    %     par.s_in_pos_mm = logpar.s_in_pos;
    % end
    % if isfield(logpar,'pos_s_pos_x') && isfield(logpar,'pos_s_pos_y')
    %     pos_s_pos_x_mm = logpar.pos_s_pos_x;
    %     pos_s_pos_y_mm = logpar.pos_s_pos_y;
    % end
    % if isempty(par.energy)
    %     energy_was_empty = 1;
    %     if isfield(logpar,'energy')
    %         par.energy = logpar.energy;
    %     end
    % else
    %     energy_was_empty = 0;
    % end
    % % if isempty(par.sample_detector_distance)
    % %     par.sample_detector_distance = logpar.sample_detector_distance;
    % % end
    % if ~exist(nexuslog_name{1},'file')
    %     % Image shape and ROI
    %     %filename = sprintf('%s%s',scan_path,ref_names{1});
    %     filename = ref_full_path{1};
    %     if ~par.raw_data
    %         %[im_raw,par.tif_info] = read_image(filename,'',[],[],[],'',par.im_trafo);
    %         [im_raw,par.tif_info] = read_image(filename,par,1);
    %         par.im_shape_raw = size(im_raw);
    %     else
    %         switch lower(cam)
    %             case 'ehd'
    %                 par.im_shape_raw = [3056 3056];
    %                 par.dtype = 'uint16';
    %             case 'kit'
    %                 par.im_shape_raw = [5120 3840];
    %                 par.dtype = 'uint16';
    %         end
    %         par.eff_pixel_size_binned = raw_bin * par.eff_pixel_size;
    %         im_raw = read_raw(filename,par.im_shape_raw,par.dtype);
    %     end
    % end
    if exist(nexuslog_name{1},'file')
        % HDF5 log
        %h5log_info = h5info(nexuslog_name{1});
        % energy,exposure time,image shape
        par.dtype = 'uint16';
        par.eff_pixel_size_binned = raw_bin * par.eff_pixel_size;
        % Image shape
        %filename = sprintf('%s%s',scan_path,ref_names{1});
        filename = ref_full_path{1};
        % mod: breaks raw data support
        % Fixed: 2019-07-10
        % CLEAN UP required
        if par.raw_data
            scanlog_name = sprintf('%s%sscan.log',scan_path,scan_name);
            if ~exist(scanlog_name,'file')
                % back up for older log file name schemes
                str = dir(sprintf('%s*scan.log',scan_path));
                scanlog_name = sprintf('%s/%s',str.folder,str.name);
            end
            [logpar,cur,cam] = p05_log(scanlog_name);
            switch lower(cam)
                case 'ehd'
                    par.im_shape_raw = [3056 3056];
                    %par.dtype = 'uint16';
                case 'kit'
                    par.im_shape_raw = [5120 3840];
                    %par.dtype = 'uint16';
            end
            %[im_raw,par.tif_info] = read_image(filename,'',[],par.tif_info,par.im_shape_raw,par.dtype,par.im_trafo);
            %else
            %[im_raw,par.tif_info] = read_image(filename,'',[],par.tif_info,[],par.dtype,par.im_trafo);
        end
        [im_raw,par.tif_info] = read_image(filename,par,1);
        par.im_shape_raw = size(im_raw);
        nexus_setup = h5info(nexuslog_name{1},'/entry/scan/setup/');
        exposure_time = h5read(nexuslog_name{1},'/entry/hardware/camera/exptime');
        magnification = h5read(nexuslog_name{1},'/entry/hardware/camera/magnification');
        pixelsize = h5read(nexuslog_name{1},'/entry/hardware/camera/pixelsize') / 1000;
        if sum(strcmpi('pos_s_stage_z',{ nexus_setup.Datasets.Name }))
            s_stage_z_setup = h5read(nexuslog_name{1},'/entry/scan/setup/pos_s_stage_z');
        end
        if sum(strcmpi('s_stage_z',{ nexus_setup.Datasets.Name }))
            s_stage_z_setup = h5read(nexuslog_name{1},'/entry/scan/setup/s_stage_z');
        end
        if sum(strcmpi('s_stage_x',{ nexus_setup.Datasets.Name }))
            s_stage_x_setup = h5read(nexuslog_name{1},'/entry/scan/setup/s_stage_x');
        end
        if sum(strcmpi('pos_s_stage_x',{ nexus_setup.Datasets.Name }))
            s_stage_x_setup = h5read(nexuslog_name{1},'/entry/scan/setup/pos_s_stage_x');
        end
        if isempty(par.eff_pixel_size)
            par.eff_pixel_size = pixelsize/magnification;
            if ~isempty(par.pixel_scaling)
                par.eff_pixel_size = par.pixel_scaling * par.eff_pixel_size;
            end
            par.eff_pixel_size_binned = raw_bin * par.eff_pixel_size;
        end
        if sum(strcmpi('pos_p05_energy',{ nexus_setup.Datasets.Name })) && energy_was_empty
            par.energy = double(h5read(nexuslog_name{1},'/entry/scan/setup/pos_p05_energy'));
        end
        if sum(strcmpi('p07_energy',{ nexus_setup.Datasets.Name })) && energy_was_empty
            par.energy = double(h5read(nexuslog_name{1},'/entry/scan/setup/p07_energy'));
        end
        if ~isempty(par.energy)
            par.energy = par.energy(end);
        end
        if isempty(par.sample_detector_distance) && sum(strcmpi('o_ccd_dist',{ nexus_setup.Datasets.Name }))
            par.sample_detector_distance = double(h5read(nexuslog_name{1},'/entry/scan/setup/o_ccd_dist')) / 1000;
        end
        if isempty(par.sample_detector_distance) && sum(strcmpi('pos_o_ccd_dist',{ nexus_setup.Datasets.Name }))
            par.sample_detector_distance = double(h5read(nexuslog_name{1},'/entry/scan/setup/pos_o_ccd_dist')) / 1000;
        end
        n_dark = h5read(nexuslog_name{1},'/entry/scan/n_dark');
        if isempty(imlogcell)
            % Get image name,key,time stamp and P3 current from log
            %[stimg_name,stimg_key,petra,petra_scan] = pp_stimg_petra(nexuslog_name,par);
            [stimg_name,stimg_key,petra,~] = pp_stimg_petra(nexuslog_name,par);

            % rotation axis
            s_rot.time = [];
            s_rot.value = [];
            for n = 1:numel(nexuslog_name)
                s_rot.time = cat(1,s_rot.time,double(h5read(nexuslog_name{n},'/entry/scan/data/s_rot/time')));
                s_rot.value = cat(1,s_rot.value,h5read(nexuslog_name{n},'/entry/scan/data/s_rot/value'));
            end
            %s_rot.time = double(h5read(nexuslog_name{1},'/entry/scan/data/s_rot/time'));
            %s_rot.value = h5read(nexuslog_name{1},'/entry/scan/data/s_rot/value');

            % lateral shift
            h5log_group = h5info(nexuslog_name{1},'/entry/scan/data/');
            if sum(strcmp('/entry/scan/data/s_stage_x',{h5log_group.Groups.Name}))
                % Read out lateral rotation axis shift form log file
                s_stage_x.time = [];
                s_stage_x.value = [];
                for n = 1:numel(nexuslog_name)
                    s_stage_x.time = cat(1,s_stage_x.time,double(h5read(nexuslog_name{n},'/entry/scan/data/s_stage_x/time')));
                    s_stage_x.value = cat(1,s_stage_x.value,h5read(nexuslog_name{n},'/entry/scan/data/s_stage_x/value'));
                end
                if ~isempty(s_stage_x.value)
                    if numel(s_stage_x.value) == numel(stimg_key.value)
                        ind = stimg_key.value == 0;
                    else
                        %ind = ~boolean(stimg_key.scan.value(logpar.n_dark+1:end));
                        ind = boolean([]);
                        for n = 1:numel(stimg_key.scan)
                            ind = cat(1,ind,~boolean(stimg_key.scan(n).value(n_dark+1:end)));
                        end
                    end
                    offset_shift_mm = s_stage_x.value(ind);
                else
                    offset_shift_mm = 0;
                end
            end
            % spiral CT translation
            h5log_group = h5info(nexuslog_name{1},'/entry/scan/data/');
            s_stage_z_str = '/entry/scan/data/s_stage_z';
        end
        % Plot PETRA current
        if par.visual_output && ~isempty(petra) && ~sum(isnan(petra.current))
            name = 'PETRA beam current from status server';
            f = figure('Name',name,'WindowState',window_state);
            x = double(petra.time - petra.time(1)) / 1000 / 60;
            y = petra.current;
            m = y ~= 0;
            plot(x(m),y(m),'.')
            xlabel('time / min')
            ylabel('current / mA')
            title(name)
            % Scan times
            if numel(stimg_key.scan.value) == numel(stimg_name.scan.time)
                t = (stimg_name.scan.time(stimg_key.scan.value == 0) - petra.time(1)) / 1000 / 60;
            else
                t = [];
            end
            if ~isempty(t)
                hold on
                y_min = min(y);
                plot(t(1),y_min,'o')
                yy = y_min + 0.05 * (max(y) - y_min);
                text(double(t(1)),yy,{'scan start','\downarrow'},'FontSize',14,'HorizontalAlignment','center')
            end
            %plot(t(end),y_min,'x')
            %text(t(end),y_min,'scan finished')
            axis tight
            drawnow
            CheckAndMakePath(fig_path)
            fig_filename = sprintf('%sfig%02u_%s.png',fig_path,f.Number,regexprep(f.Name,'\ |:','_')) ;
            saveas(f,fig_filename);
        end

        %% Lateral shift
        if numel(offset_shift_mm) > 1
            offset_shift_mm = offset_shift_mm(par.proj_range);
        end
        if numel(offset_shift_mm) && abs(std(offset_shift_mm)) * 1000 > 1
            % Shift or static position            
            if std(offset_shift_mm) > 10 * eps
                offset_shift = 1e-3 / par.eff_pixel_size * offset_shift_mm;
                par.s_in_pos = 1e-3 / par.eff_pixel_size * par.s_in_pos_mm;
                offset_shift_min = min(offset_shift(:)) ;
                offset_shift = 1 + offset_shift - offset_shift_min;
                % Overwrite lateral shift if offset shift is provided as parameter
                if isequal(std(offset_shift),0) && ~isempty(tomo.rot_axis_offset_shift) %&& isscalar(tomo.rot_axis_offset_shift)
                    %offset_shift = tomo.rot_axis_offset_shift / raw_bin * (0:par.num_proj_used) / par.num_proj_used;
                    offset_shift = tomo.rot_axis_offset_shift / raw_bin;% * (0:par.num_proj_used) / par.num_proj_used;
                end
                % Tranform to integer pixel-wise shifts
                tmp = offset_shift;
                offset_shift = round(offset_shift);
                os_std = std(offset_shift - tmp);
                if os_std > 1e-2
                    cprintf('Red','\nOffset shift not on integer pixel scale: %f (%f %f)',os_std,offset_shift([1,2])-tmp([1,2]))
                end
                % Lateral scanning
                % position index extracted by jump in offset_shift
                scan_position_index = zeros(size(offset_shift));
                num_scan_pos = 1;
                if offset_shift_mm(1) < max(offset_shift_mm)
                    scan_dir = 1;
                else
                    scan_dir = -1;
                end
                scan_position_index(1) = num_scan_pos;
                for nn = 2:numel(offset_shift)
                    if abs(offset_shift(nn) - offset_shift(nn-1)) > 601
                        num_scan_pos = num_scan_pos + 1;
                    end
                    scan_position_index(nn) = num_scan_pos;
                end
                % Relative scan position without lateral offset
                scan_position = zeros(size(offset_shift));
                for nn = 1:num_scan_pos
                    m = scan_position_index == nn;
                    scan_position(m) = min(offset_shift(m)) - 1;
                end
                offset_shift = offset_shift - scan_position;
                % Scale position because of binning for tomo reco
                scan_position = scan_position - mean(scan_position);
                scan_position = 1 / raw_bin * scan_position;
                if tomo.flip_scan_position
                    scan_position = - scan_position;
                end
                fprintf(' \n scan positions inclusive multiple turns : %u',num_scan_pos)
                num_lateral_pos = unique(scan_position);
                num_lateral_pos = numel(num_lateral_pos);
                fprintf(' \n lateral scan positions : %u',num_lateral_pos)
                if isempty(tomo.vol_size)
                    tomo.vol_size = [num_lateral_pos num_lateral_pos num_lateral_pos num_lateral_pos 1 1] .* [-0.5 0.5 -0.5 0.5 -0.5 0.5];
                    fprintf('\nSetting tomo.vol_size = [')
                    fprintf(' %g',tomo.vol_size)
                    fprintf(']')
                    % infer (relative) vol_shape from vol_size if vol_shape is empty,but vol_size is given
                    if isempty(tomo.vol_shape)
                        tomo.vol_shape = tomo.vol_size(2:2:end) - tomo.vol_size(1:2:end);
                        fprintf('\nSetting tomo.vol_shape = [')
                        fprintf(' %g',tomo.vol_shape)
                        fprintf(']')
                    end
                end
                % Print info
                x0 = min(offset_shift_mm);
                x1 = max(offset_shift_mm);
                dx = x1 - x0;
                fprintf(' \n absolute lateral shift / micron : [min max dx] = [%f %f %f] ',x0,x1,dx)
                x0 = min(offset_shift);
                x1 = max(offset_shift);
                dx = x1 - x0;
                fprintf(' \n relative lateral shift / pixel (shaking): [min max dx] = [%f %f %f] ',x0,x1,dx)
                fprintf(' \n scan direction : %.0f',scan_dir)
                if scan_dir == 1
                    fprintf(' (from right to left w.r.t. to the sample)')
                elseif scan_dir == -1
                    fprintf(' (from left to right w.r.t. to the sample)')
                end
                % Plot offset shift%
                if par.visual_output %&& std(scan_position) ~= 0
                    f = figure('Name','rotation axis offset shift','WindowState',window_state);

                    yyaxis left
                    plot(offset_shift_mm,'.')
                    title(sprintf('rotation axis offset shift. effective pixel size binned: %.2f micron,scan dir: %u',par.eff_pixel_size_binned * 1e6,scan_dir))

                    axis tight
                    xlabel('projection number')
                    ylabel('absolute lateral shift (s stage x) / mm')

                    yyaxis right
                    plot(scan_position,'.')
                    ylabel('relative scan position / pixel')
                    ymin = min(scan_position) - dx / 2 / raw_bin;
                    ymax = max(scan_position) + dx / 2 / raw_bin;
                    ylim([ymin ymax])

                    legend({ 's stage x','relative scan position' })
                    drawnow
                    CheckAndMakePath(fig_path)
                    fig_filename = sprintf('%sfig%02u_%s.png',fig_path,f.Number,regexprep(f.Name,'\ |:','_'));
                    saveas(f,fig_filename);
                end
            end % if std(offset_shift_mm)
        else
            par.crop_proj = 0;
        end % if numel(s_stage_x.value)
        if isempty(par.sample_detector_distance)
            par.sample_detector_distance = logpar.sample_detector_distance;
        end
        %% Vertical shift
        if ~isempty(s_stage_z_str) && sum(strcmp(s_stage_z_str,{h5log_group.Groups.Name}))
            % Read out vertical shift form HDF5 log file
            s_stage_z.time = [];
            s_stage_z.value = [];
            for n = 1:numel(nexuslog_name)
                s_stage_z.time = cat(1,s_stage_z.time,double(h5read(nexuslog_name{n},[s_stage_z_str '/time'])));
                s_stage_z.value = cat(1,s_stage_z.value,h5read(nexuslog_name{n},[s_stage_z_str '/value']));
            end

            % Static or shift?
            if numel(s_stage_z.value) % is not empty
                switch numel(s_stage_z.value)
                    case par.num_proj_found
                        vert_shift_micron = s_stage_z.value;
                    case par.num_proj_found + par.num_dark + par.num_ref_found
                        m = stimg_key.scan.value == 0 ;
                        vert_shift_micron = s_stage_z.value(m);
                    case par.num_proj_found + par.num_ref_found
                        %m = stimg_key.value(logpar.n_dark + 1:end) == 0 ;
                        %m = stimg_key.value(par.num_proj_found + 1:end) == 0 ;
                        vert_shift_micron = s_stage_z.value(round(par.num_proj_found/2) + (1:par.num_proj_found));
                    otherwise
                        vert_shift_micron = s_stage_z.value(1:par.num_proj_found);
                end
                vert_shift_micron = vert_shift_micron(par.proj_range);
                if abs(std(SubtractMean(vert_shift_micron))) > 1e-3
                    %vert_shift_micron = vert_shift_micron(par.proj_range);
                    % Check
                    vert_shift = vert_shift_micron * 1e-3 / par.eff_pixel_size_binned;
                    vert_shift = SubtractMean(vert_shift);
                    fprintf('\n vertical shift absolute / micron : [%g %g]',min(vert_shift_micron),max(vert_shift_micron))
                    fprintf('\n vertical shift relative / binned pixel : [%g %g] ',min(vert_shift),max(vert_shift))
                    dz_micron = max(vert_shift_micron) - min(vert_shift_micron);
                    dz = max(vert_shift) - min(vert_shift);
                    fprintf('\n vertical shift diff : %g micron,%g binned pixel',dz_micron,dz)
                    fprintf('\n vertical shift / #proj / unbinned pixel : %g',dz / par.num_proj_found * raw_bin);
                    % Plot vertical shift
                    if par.visual_output && std(vert_shift) ~= 0
                        name = 'spiral scan: vertical shift';
                        f = figure('Name',name,'WindowState',window_state);
                        plot(vert_shift,'.')
                        title(name)
                        axis equal tight
                        xlabel('projection number')
                        ylabel('vertical / pixel')
                        legend(sprintf('shift / #proj / pixel %f',dz / par.num_proj_found))
                        drawnow
                        CheckAndMakePath(fig_path)
                        fig_filename = sprintf('%sfig%02u_%s.png',fig_path,f.Number,regexprep(f.Name,'\ |:','_'));
                        saveas(f,fig_filename);
                    end
                end % if std(vert_shift_micron)
            end % if numel(s_stage_z.value)
        end %if sum(strcmp('/entry/scan/data/s_stage_z',{a.Groups.Name}))
        %% Ring current
        if ~isempty(petra) && numel(stimg_key.scan.value) == numel(stimg_name.scan.time)
            X = double(petra.time);
            V = double(petra.current);
            Xq = double(stimg_name.scan.time);
            extrap_val = median(V);
            if isequal(X,Xq)
                stimg_name.current = V;
            else
                stimg_name.current = (interp1(X,V,Xq,'next',extrap_val) + interp1(X,V,Xq + double(exposure_time),'previous',extrap_val)) / 2;
            end
            cur_ref_val = stimg_name.current(stimg_key.scan.value == 1);
            cur_ref_name = stimg_name.scan.value(stimg_key.scan.value == 1);
            cur_ref_time = stimg_name.scan.time(stimg_key.scan.value == 1);
            re = regexp(cur_ref_name{1},'\d{7,7}');
            len = 7;
            if isempty(re)
                re = regexp(cur_ref_name{1},'\d{6,6}');
                len = 6;
            end
            if isempty(re)
                re = regexp(cur_ref_name{1},'\d{5,5}');
                len = 5;
            end
            if isempty(re)
                re = regexp(cur_ref_name{1},'\d{4,4}');
                len = 4;
            end
            if isempty(re)
                len = [];
            end
            if numel(re) >= 1
                re = re(end);
                imtype_str_flag = re;
            elseif strcmpi(ref_names{1}(end-6:end-4),'ref')
                imtype_str_flag = '64'; % 0
            elseif strcmpi(ref_names{1}(end-11:end-9),'ref')
                imtype_str_flag = '119'; % 1
            end
            for nn = numel(cur_ref_name):-1:1
                cur.ref(nn).val = cur_ref_val(nn);
                cur.ref(nn).name = cur_ref_name{nn};
                cur.ref(nn).time = cur_ref_time(nn);
                if imtype_str_flag >= 0
                    if ~isempty(len)
                        cur.ref(nn).ind = str2double(cur.ref(nn).name(imtype_str_flag + (0:len-1)));
                    else
                        cur.ref(nn).ind = str2double(cur.ref(nn).name(imtype_str_flag + (0:5)));
                    end
                elseif imtype_str_flag == -1
                    cur.ref(nn).ind = str2double(cur.ref(nn).name(end-12:end-8));
                elseif imtype_str_flag == -2
                    cur.ref(nn).ind = str2double(cur.ref(nn).name(end-7:end-4));
                end
            end
            cur_proj_val = stimg_name.current(stimg_key.scan.value == 0);
            cur_proj_name = stimg_name.scan.value(stimg_key.scan.value == 0);
            cur_proj_time = stimg_name.scan.time(stimg_key.scan.value == 0);
            for nn = numel(cur_proj_name):-1:1
                cur.proj(nn).val = cur_proj_val(nn);
                cur.proj(nn).name = cur_proj_name{nn};
                cur.proj(nn).time = cur_proj_time(nn);
                if imtype_str_flag >= 0
                    if ~isempty(len)
                        cur.proj(nn).ind = str2double(cur.proj(nn).name(imtype_str_flag + (0:len-1)));
                    else
                        cur.proj(nn).ind = str2double(cur.proj(nn).name(imtype_str_flag + (0:5)));
                    end
                elseif imtype_str_flag == -1
                    cur.proj(nn).ind = str2double(cur.proj(nn).name(end-12:end-8));
                elseif imtype_str_flag == -2
                    cur.proj(nn).ind = str2double(cur.proj(nn).name(end-7:end-4));
                end
            end
        end
        %t0 = min([cur_proj_time; cur_ref_time]);
    else % if ~exist(nexuslog_name,'file')
        scanlog_name = sprintf('%s%sscan.log',scan_path,scan_name);
        if ~exist(scanlog_name,'file')
            % back up for older log file name schemes
            str = dir(sprintf('%s*scan.log',scan_path));
            scanlog_name = sprintf('%s/%s',str.folder,str.name);
        end
        [logpar,cur,cam] = p05_log(scanlog_name);
        if isempty(par.eff_pixel_size)
            par.eff_pixel_size = logpar.eff_pixel_size;
            if ~isempty(par.pixel_scaling)
                par.eff_pixel_size = par.pixel_scaling * par.eff_pixel_size;
            end
            par.eff_pixel_size_binned = raw_bin * par.eff_pixel_size;
        end
        if isempty(exposure_time)
            exposure_time = logpar.exposure_time;
        end
        if isfield(logpar,'s_in_pos')
            par.s_in_pos_mm = logpar.s_in_pos;
        end
        if isfield(logpar,'pos_s_pos_x') && isfield(logpar,'pos_s_pos_y')
            pos_s_pos_x_mm = logpar.pos_s_pos_x;
            pos_s_pos_y_mm = logpar.pos_s_pos_y;
        end
        if isempty(par.energy)
            if isfield(logpar,'energy')
                par.energy = logpar.energy;
            end
        end
        % if isempty(par.sample_detector_distance)
        %     par.sample_detector_distance = logpar.sample_detector_distance;
        % end
        if ~exist(nexuslog_name{1},'file')
            % Image shape and ROI
            %filename = sprintf('%s%s',scan_path,ref_names{1});
            filename = ref_full_path{1};
            if ~par.raw_data
                %[im_raw,par.tif_info] = imagemage(filename,'',[],[],[],'',par.im_trafo);
                [im_raw,par.tif_info] = read_image(filename,par,1);
                par.im_shape_raw = size(im_raw);
            else
                switch lower(cam)
                    case 'ehd'
                        par.im_shape_raw = [3056 3056];
                        par.dtype = 'uint16';
                    case 'kit'
                        par.im_shape_raw = [5120 3840];
                        par.dtype = 'uint16';
                end
                par.eff_pixel_size_binned = raw_bin * par.eff_pixel_size;
                im_raw = read_raw(filename,par.im_shape_raw,par.dtype);
            end
        end
        if par.crop_proj
            if ~isempty(tomo.rot_axis_offset_shift)
                offset_shift = tomo.rot_axis_offset_shift;
                %offset_shift = SubtractMean(offset_shift);
                offset_shift = 1e-3 / par.eff_pixel_size * offset_shift;
                offset_shift_min = min(offset_shift(:)) ;
                offset_shift = 1 + offset_shift - offset_shift_min;
                % Tranform to integer pixel-wise shifts
                tmp = offset_shift;
                offset_shift = round(offset_shift);
                os_std = std(offset_shift - tmp);
                if os_std > 1e-2
                    cprintf('Red',sprintf('\nOffset shift not on integer pixel scale: %f\n',os_std))
                end
            end
        end
    end % if ~exist(nexuslog_name,'file')
    %% Raw ROI
    par = set_raw_roi(par,im_raw,scan_path,fig_path,ref_full_path,dark_names);
    %% Print info
    %im_roi = read_image(filename,'',par.raw_roi,par.tif_info,par.im_shape_raw,par.dtype,par.im_trafo);
    par.im_roi = read_image(filename,par);
    im_shape_roi = size(par.im_roi);
    im_shape_binned1 = floor(size(par.im_roi,1) / raw_bin);
    im_shape_binned2 = floor(size(par.im_roi,2) / raw_bin);
    par.im_shape_binned1 = im_shape_binned1;
    par.im_shape_binned2 = im_shape_binned2;
    fprintf('\n energy : %.1f keV',par.energy / 1e3)
    par.wave_length = E_to_lambda(par.energy);
    fprintf('\n wave length lambda: %.2g angstrom',par.wave_length * 10^10)
    fprintf('\n wave number 2*pi/lambda %f 1/nm:',2 * pi / par.wave_length / 10^9)
    fprintf('\n distance sample dector : %.1f mm',par.sample_detector_distance * 1000)
    fprintf('\n effective pixel size unbinned : %.2f micron', par.eff_pixel_size * 1e6)
    fprintf('\n effective pixel size binned: %.2f micron', par.eff_pixel_size_binned * 1e6)
    par.fresnel_number = par.eff_pixel_size_binned^2 / (par.wave_length * par.sample_detector_distance);
    fprintf('\n Fresnel number =  pixelsize^2 / lambda / z: %.2f', par.fresnel_number)
    par.exposure_time = exposure_time;
    pf = 2 * pi /par.fresnel_number;
    zercro = pf/pi/4;
    fprintf('\n Fresnel prefactor = 2*pi*lambda*z/pixelsize_m^2 = %06.3f',pf)
    fprintf('\n maximum argument of the sine at [|xi| |eta|] = [1/2 1/2]: %5.3f',pf/4)
    fprintf('\n zero crossings - 1 (without the central one at xi = eta = 0): %5.3f',zercro);
    fprintf('\n exposure time : %g ms',exposure_time);
    fprintf('\n exposure time projections : %g s',exposure_time * par.num_proj_found / 1000);
    fprintf('\n exposure time flats fields : %g s',exposure_time * par.num_ref_found / 1000);
    exp_total = exposure_time * (par.num_proj_found + par.num_ref_found + par.num_dark);
    fprintf('\n exposure time total : %g s = %g min',exp_total / 1000,exp_total / 1000 / 60);
    fprintf('\n image size : %.2f mm x %.2f mm',par.im_shape_raw * par.eff_pixel_size *1e3)
    fprintf('\n image shape : %u x %u = %.1g pixels',par.im_shape_raw,prod(par.im_shape_raw))
    fprintf('\n image shape roi : %u x %u = %.1g pixels',im_shape_roi,numel(par.im_roi))
    numel_im_roi_binned = im_shape_binned1 * im_shape_binned2;
    fprintf('\n image shape roi binned : %u x %u = %.1g pixels',im_shape_binned1,im_shape_binned2,numel_im_roi_binned)
    fprintf('\n raw binning factor : %u',raw_bin)
    if par.use_gpu_in_parfor
        % Limit by GPU pixel filter: 2 * single + 2 * uint16 + 1 * logical
        gpu_mem_requ_per_im = prod(im_shape_roi + 2 * pixel_filter_radius) * (4 + 4 + 2 + 2 + 1);
        poolsize_max_gpu = floor(par.poolsize_gpu_limit_factor * min(mem_avail_gpu) / gpu_mem_requ_per_im);
        %% TODO: Check max poolsize
        %poolsize_max_gpu = max([poolsize_max_gpu,numel(par.gpu_index)]);
        poolsize_max_gpu = poolsize_max_gpu * numel(par.gpu_index);
        fprintf(' \n estimated GPU memory required per image for pixel filtering : %g MiB',gpu_mem_requ_per_im / 1024^2)
        fprintf(' \n GPU poolsize limit factor : %g',par.poolsize_gpu_limit_factor)
        fprintf(' \n GPU memory induced maximum poolsize : %u',poolsize_max_gpu)
    else
        poolsize_max_gpu = par.poolsize;
    end
    %% Dark field %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t = toc;
    fprintf('\nProcessing %u dark fields.',par.num_dark)
    dark = zeros([im_shape_roi,par.num_dark],'uint16');
    fprintf('\n allocated memory: %.2f MiB',Bytes(dark,2))
    filt_pix_par.threshold_hot = pixel_filter_threshold_dark(1);
    filt_pix_par.threshold_dark = pixel_filter_threshold_dark(2);
    filt_pix_par.medfilt_neighboorhood = pixel_filter_radius;
    filt_pix_par.filter_dead_pixel = 1;
    filt_pix_par.filter_Inf = 0;
    filt_pix_par.filter_NaN = 0;
    filt_pix_par.verbose = 0;
    filt_pix_par.use_gpu = par.use_gpu_in_parfor;
    read_image_par = @(filename) read_image(filename,par);
    FilterPixel_par = @(im_int) FilterPixel(im_int,filt_pix_par);
    % par loop dark
    parfor (nn = 1:par.num_dark,poolsize_max_gpu)
        % Read image
        %filename = sprintf('%s%s',scan_path,dark_names{nn});
        filename = dark_full_path{nn};
        im_int = read_image_par(filename);
        % Remove large outliers,hack for KIT camera chip artefacts. Assume
        % Gaussian distribution and set all value above mean + 4 * std
        % (99.994 of values lie within 4 std).
        im_float = single(im_int(:));
        im_mean = mean(im_float);
        im_std = std(im_float);
        im_int(im_int > im_mean + 4*im_std) = uint16(im_mean);
        % Filter pixels
        im_int = FilterPixel_par(im_int);
        % Assign image to stack
        dark(:,:,nn) = im_int;
    end
    stats.darks_min = min(dark(:));
    stats.darks_max = max(dark(:));
    % Reject dark images which are all zero
    darks_to_use = zeros(1,par.num_dark,'logical');
    parfor nn = 1:par.num_dark
        darks_to_use(nn) = boolean(max2(dark(:,:,nn)) );
    end
    % Median/Mean dark
    dark_median = squeeze(median(dark(:,:,darks_to_use),3));
    %dark_mean = squeeze(mean(dark(:,:,darks_to_use),3));
    stat.dark_med_min = min(dark(:));
    stat.dark_med_max = max(dark(:));
    dark = dark_median;
    % Binned dark
    dark_median_binned = 1 / raw_bin^2 * Binning(dark_median,raw_bin);
    %dark_mean_binned = 1 / raw_bin^2 * Binning(dark_mean,raw_bin);
    fprintf('\n duration : %.1f s',toc-t)
    fprintf('\n min/max of all darks : %g %g',stats.darks_min,stats.darks_max);
    fprintf('\n min/max of median dark : %g %g',stat.dark_med_min,stat.dark_med_max);
    % Save dark
    CheckAndMakePath(im_path3)
    write32bitTIFfromSingle(sprintf('%sdark_median_binned.tif',im_path3),rot90(dark_median_binned))
    %write32bitTIFfromSingle(sprintf('%sdark_mean_binned.tif',im_path3),rot90(dark_mean_binned))
    % Fig: raw + dark field
    if par.visual_output
        h1 = figure('Name','data and flat-and-dark-field correction','WindowState',window_state);
        subplot(2,3,1)
        imsc1(dark_median_binned);
        title(sprintf('median dark field'))
        colorbar
        axis equal tight
        xticks('auto'),yticks('auto')
        drawnow
    end
    %% Flat field %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t = toc;
    fprintf('\nProcessing %u flat fields.',par.num_ref_used)
    % Correlation roi area
    flat_corr_area1 = IndexParameterToRange(image_correlation.area_width,im_shape_roi(1));
    flat_corr_area2 = IndexParameterToRange(image_correlation.area_height,im_shape_roi(2));
    % Preallocation
    flat = zeros([im_shape_binned1,im_shape_binned2,par.num_ref_used],'single');
    flat_min = zeros([1 par.num_ref_used],'single');
    flat_max = zeros([1 par.num_ref_used],'single');
    flat_mean = zeros([1 par.num_ref_used],'single');
    flat_std = zeros([1 par.num_ref_used],'single');
    num_zeros = zeros(1,par.num_ref_used,'single');
    flat_corr_area_width_binned = floor(numel(flat_corr_area1) / 2 / raw_bin);
    flat_corr_area_height_binned = floor(numel(flat_corr_area2) / 2 / raw_bin);
    roi_flat = zeros(flat_corr_area_width_binned,flat_corr_area_height_binned,par.num_ref_used,'single');
    fprintf('\n allocated memory: %.2f GiB',Bytes(flat,3))
    refs_to_use = zeros(1,size(flat,3),'logical');
    filt_pix_par.threshold_hot = pixel_filter_threshold_flat(1);
    filt_pix_par.threshold_dark = pixel_filter_threshold_flat(2);
    % image correlation filter function for parloop
    switch image_correlation.filter_type
        case 'median'
            imcf_filter = @(im) medfilt2(im,image_correlation.filter_parameter{:});
        case 'wiener'
            imcf_filter = @(im) wiener2(im,image_correlation.filter_parameter{:});
        otherwise
            imcf_filter = @(im) im;
    end
    % Parallel loop over refs
    startS = ticBytes(gcp);
    % read_image_par = read_image(filename,'',par.raw_roi,par.tif_info,par.im_shape_raw,par.dtype,par.im_trafo);
    %ref_names_mat = NameCellToMat(ref_names(par.ref_range));
    read_image_par = @(filename) read_image(filename,par);
    % par loop refs
    parfor (nn = 1:par.num_ref_used,poolsize_max_gpu)
        % Read
        filename = ref_full_path{nn};
        im_int = read_image_par(filename);
        % Filter pixel
        im_int = FilterPixel(im_int,filt_pix_par);
        % Dark field correction
        im_int = im_int - dark;
        % Correlation ROI
        roi_flat(:,:,nn) = 1 / 2 / raw_bin^2 * Binning(im_int(flat_corr_area1,flat_corr_area2),2 * raw_bin);
        % Correlation ROI filter
        roi_flat(:,:,nn) = imcf_filter(roi_flat(:,:,nn)); %#ok<*PFBNS>
        % Binning
        im_float_binned = Binning(im_int,raw_bin) / raw_bin^2;
        % Count for zeros
        num_zeros(nn) =  nnz(im_float_binned(:) < 1 );
        % Discard if any pixel is zero.
        refs_to_use(nn) = ~boolean(num_zeros(nn) );
        % Assign image to stack
        flat(:,:,nn) = im_float_binned;
        % Statistics
        flat_min(nn) = min2(im_float_binned);
        flat_max(nn) = max2(im_float_binned);
        flat_mean(nn) = mean2(im_float_binned);
        flat_std(nn) = std2(im_float_binned);
    end
    if ~sum(refs_to_use)
        cprintf('Red','\nWARNING: All flat fields discarded\n')
        dbstop
    end
    toc_bytes.read_flat = tocBytes(gcp,startS);
    % Plot image statistics
    if par.visual_output
        name = 'image statistics: flat fields';
        his = figure('Name',name,'WindowState',window_state);
        % flat min
        subplot(2,2,1);
        plot(flat_min,'.')
        axis tight
        xlabel('image no.')
        title('min')
        % flat max
        subplot(2,2,2);
        plot(flat_max,'.')
        axis tight
        xlabel('image no.')
        title('max')
        % flat mean
        subplot(2,2,3);
        plot(flat_mean,'.')
        axis tight
        xlabel('image no.')
        title('mean')
        % flat std
        subplot(2,2,4);
        plot(flat_std,'.')
        axis tight
        xlabel('image no.')
        title('std')
        drawnow
        CheckAndMakePath(fig_path)
        fig_filename = sprintf('%sfig%02u_%s.png',fig_path,his.Number,regexprep(his.Name,'\ |:','_'));
        saveas(his,fig_filename);
    end
    % Delete empty refs
    zz = ~refs_to_use;
    if sum(zz(:))
        fprintf(' Deleting empty flats')
        flat(:,:,~refs_to_use) = [];
    end
    % min/max values before dark field subtraction and ring current normalization
    stats.flat_min = min(flat(:));
    stats.flat_max = max(flat(:));
    % Ring current normalization
    if par.ring_current_normalization && isempty(par.ref_path)
        ref_ind_from_filenames = CellString2Vec(ref_names(par.ref_range));
        ref_ind_from_log = [cur.ref(par.ref_range).ind];
        if isequal(ref_ind_from_filenames,ref_ind_from_log)
            ref_rc = [cur.ref(par.ref_range).val];
            %ref_t = ([cur.ref(par.ref_range).time] - t0) / 1000 / 60;
            ref_ind = [cur.ref(par.ref_range).ind];
            ref_rcm = mean(ref_rc(:));
            scale_factor = 100 ./ shiftdim(ref_rc(refs_to_use),-1);
            flat = fun_times(flat,scale_factor);
            if par.visual_output
                hrc = figure('Name','PETRA III beam current: Interpolation at image time stamps','WindowState',window_state);
                subplot(1,1,1);
                %plot(ref_rc(:),'.')
                %plot(ref_t(:),ref_rc(:),'.')
                plot(ref_ind(:),ref_rc(:),'.')
                axis tight
                title(sprintf('ring current: flat fields'))
                legend(sprintf('mean: %.2f mA',ref_rcm))
                drawnow
            end
        else
            cprintf('Red','\n Flat fields not normalized by ring current. Names read from directory and log-file are inconsistent.')
        end
    else
        cprintf('Red','\n Flat fields not normalized by ring current.')
    end
    stats.flat_corrected_min = min(flat(:));
    stats.flat_corrected_max = max(flat(:));
    nn =  nnz(flat(:) < 1);
    if nn > 0
        cprintf('Red','\nWARNING: Flat field contains %u zeros\n',nn)
    end
    nn = sum(~refs_to_use(:));
    par.num_ref_used = par.num_ref_used - nn;
    fprintf('\n duration : %.1f s',toc-t)
    fprintf('\n min/max of all flats : %6g %6g',stats.flat_min,stats.flat_max);
    fprintf('\n min/max of all corrected flats : %6g %6g',stats.flat_corrected_min,stats.flat_corrected_max);
    PrintVerbose(nn,'\n discarded empty refs : %u,%.2f%%',nn,100 * nn / par.num_ref_found)
    if sum(num_zeros)
        fprintf('\n flat fields with zeros :')
        % print #zeros if not all pixels are zero
        for nn = 1:numel(num_zeros)
            if num_zeros(nn) ~= 0
                if isequal(num_zeros(nn),numel_im_roi_binned)
                    fprintf(' %u',nn)
                else
                    fprintf(' %u:%u',nn,num_zeros(nn))
                end
            end
        end
    end
    % Save flat images
    num_flat12 = round(size(flat,3) / 2) ;
    write32bitTIFfromSingle(sprintf('%sflat_dark_subtracted_beamcurrent_corrected_binned_%06u.tif',im_path3,1),rot90(flat(:,:,1)))
    write32bitTIFfromSingle(sprintf('%sflat_dark_subtracted_beamcurrent_corrected_binned_%06u.tif',im_path3,num_flat12),rot90(flat(:,:,num_flat12))) ;
    write32bitTIFfromSingle(sprintf('%sflat_dark_subtracted_beamcurrent_corrected_binned_%06u.tif',im_path3,size(flat,3)),rot90(flat(:,:,end)))
    write32bitTIFfromSingle(sprintf('%sflat_dark_subtracted_beamcurrent_corrected_binned_mean3.tif',im_path3),rot90(mean(flat,3)))
    %% Figure: Flat field
    if par.visual_output
        % Show flat field
        if exist('h1','var') && isvalid(h1)
            figure(h1)
        else
            h1 = figure('Name','data and flat-and-dark-field correction','WindowState',window_state);
        end
        subplot(2,3,2)
        imsc1(flat(:,:,1))
        title(sprintf('flat field #1'))
        colorbar
        axis equal tight
        drawnow
        % Correlation area
        if ~sum(strcmp(image_correlation.method,{'none','median','mean'}))
            h_corr_roi = figure('Name','image correlation roi','WindowState',window_state);
            % mean flat
            subplot(2,2,1)
            im = mean(flat,3);
            imsc1(im)
            hold on
            x1 = ceil(flat_corr_area1(1) / raw_bin);
            y1 = im_shape_binned2 - floor(flat_corr_area2(end) / raw_bin);
            rectangle('position',[ x1 y1 2 * flat_corr_area_width_binned 2 * flat_corr_area_height_binned],'EdgeColor','r')
            title(sprintf('mean flat field\n(vertical coordiante values flipped)'))
            colorbar
            axis equal tight
            % ROI flat
            subplot(2,2,2)
            imsc1(roi_flat(:,:,1))
            title(sprintf('flat field ROI'));
            ylabel(sprintf('image correlation area:\nunbinned (relative) coordinates\nwidth = %g %g\nheight = %g %g',image_correlation.area_width,image_correlation.area_height))
            colorbar
            axis equal tight
            drawnow
        end
    end
    %% Projections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t = toc;
    fprintf('\nProcessing %u projections.',par.num_proj_used)
    %img_names_mat = NameCellToMat(proj_names(par.proj_range));
    filt_pix_par.threshold_hot = pixel_filter_threshold_proj(1);
    filt_pix_par.threshold_dark = pixel_filter_threshold_proj(2);
    % Display first raw image
    if par.visual_output
        if exist('h1','var') && isvalid(h1)
            figure(h1)
        else
            h1 = figure('Name','data and flat-and-dark-field correction','WindowState',window_state);
        end
        filename = proj_full_path{1};
        im_int = read_image(filename,par);
        im_int = FilterPixel(im_int,filt_pix_par);
        raw1 = Binning(im_int,raw_bin) / raw_bin^2;
        subplot(2,3,3)
        imsc1(raw1)
        title(sprintf('raw proj #1'))
        colorbar
        axis equal tight
        drawnow
    end
    % Get absolut filter thresholds from percentage-wise pixel filtering
    % of 1st,middle,and last projection to speed up processing
    if filt_pix_par.threshold_hot  < 1 || filt_pix_par.threshold_dark < 0.5
        %filename = sprintf('%s%s',scan_path,img_names_mat(par.num_proj_used,:));
        filename = proj_full_path{par.num_proj_used};
        im_int = read_image(filename,par);
        [~,ht(3),dt(3)] = FilterPixel(im_int,filt_pix_par);
        %filename = sprintf('%s%s',scan_path,img_names_mat(1,:));
        filename = proj_full_path{1};
        im_int = read_image(filename,par);
        [~,ht(2),dt(2)] = FilterPixel(im_int,filt_pix_par);
        %filename = sprintf('%s%s',scan_path,img_names_mat(round(par.num_proj_used/2),:));
        filename = proj_full_path{round(par.num_proj_used/2)};
        im_int = read_image(filename,par);
        [~,ht(1),dt(1)] = FilterPixel(im_int,filt_pix_par);
        filt_pix_par.threshold_hot  = median(ht);
        filt_pix_par.threshold_dark = median(dt);
    end
    % Preallocation
    roi_proj = zeros(floor(numel(flat_corr_area1) / 2 / raw_bin),floor(numel(flat_corr_area2) / 2 / raw_bin),par.num_proj_used,'single');
    % Lateral shift indices
    if ~par.crop_proj || isscalar(offset_shift)
        x0 = ones(1,par.num_proj_used);
        x1 = im_shape_roi(1) * x0;
    else
        x0 = offset_shift';
        x1 = x0 + im_shape_roi(1) - max(x0);
    end
    % Parallel loop over projections
    %read_image_par = @(filename) read_image(filename,'',par.raw_roi,par.tif_info,par.im_shape_raw,par.dtype,par.im_trafo);
    read_image_par = @(filename) read_image(filename,par);
    startS = ticBytes(gcp);
    % Preallocation
    im_shape_cropbin1 = floor((x1(1) - x0(1) + 1) / raw_bin);
    par.im_shape_cropbin1 = im_shape_cropbin1;
    proj = zeros(im_shape_cropbin1,im_shape_binned2,par.num_proj_used,'single');
    proj_min = zeros([1 par.num_proj_used],'single');
    proj_max = zeros([1 par.num_proj_used],'single');
    proj_mean = zeros([1 par.num_proj_used],'single');
    proj_std = zeros([1 par.num_proj_used],'single');
    %num_zeros = zeros(1,par.num_proj_used,'single');
    fprintf('\n allocated memory: %.2f GiB',Bytes(proj,3))
    projs_to_use = zeros(1,size(proj,3),'logical');
    dep = par.delete_empty_projections;
    dofc = 0;
    flat_m = 0;
    switch image_correlation.method
        case 'mean'
            flat_m = mean(flat,3);
            dofc = 1;
        case 'median'
            flat_m = median(flat,3);
            dofc = 1;
    end
    % par loop projections
    parfor (nn = 1:par.num_proj_used,poolsize_max_gpu)
        %im = proj(:,:,nn);
        % Read projection
        %filename = sprintf('%s%s',scan_path,img_names_mat(nn,:));
        filename = proj_full_path{nn};
        %im_int = read_image(filename,'',par.raw_roi,par.tif_info,par.im_shape_raw,par.dtype,par.im_trafo);
        im_int = read_image_par(filename);
        % Filter pixel
        im_int = FilterPixel(im_int,filt_pix_par);
        % Dark field correction
        im_int = im_int - dark;
        % Correlation ROI
        roi_proj(:,:,nn) = 1 / 2 / raw_bin^2 * Binning(im_int(flat_corr_area1,flat_corr_area2),2 * raw_bin);
        % Correlation ROI filter
        roi_proj(:,:,nn) = imcf_filter(roi_proj(:,:,nn));
        % Remove lateral shift & Binning
        xx = x0(nn):x1(nn);
        im_float_binned = Binning(im_int(xx,:),raw_bin) / raw_bin^2;
        % Count zeros
        %num_zeros(nn) = nnz(im_float_binned(:) < 1);
        num_zeros = nnz(im_float_binned(:) < 1);
        % Reject image if any pixel is zero
        if dep
            %projs_to_use(nn) = ~boolean(num_zeros(nn));
            projs_to_use(nn) = ~boolean(num_zeros);
        else
            projs_to_use(nn) = 1;
        end
        % % Binned shift (shift,not first pixel)
        % shift = (x0(nn) - 1) / raw_bin;
        % shift_int = floor(shift);
        % shift_sub = shift - shift_int;
        %
        % % shift flat
        % if mod(shift_sub,1) ~= 0
        %     % crop flat at integer shift,then shift subpixel
        %     xx = shift_int + (1:im_shape_cropbin1+1);
        %     flat_median_shifted = imtranslate(flat_m(xx,:),[0 -shift_sub],'linear');
        % else
        %     xx = shift_int + (1:im_shape_cropbin1);
        %     flat_median_shifted = flat_m(xx,:);
        % end
        %
        % % flat field correction
        % p = proj(:,:,nn);
        % p = p ./ flat_median_shifted(1:im_shape_cropbin1,:) ;
        % % Reassign
        if dofc
            im_float_binned = im_float_binned ./ flat_m;
        end
        % Assign image to stack
        proj(:,:,nn) = im_float_binned;
        % Statistics
        proj_min(nn) = min2(im_float_binned);
        proj_max(nn) = max2(im_float_binned);
        proj_mean(nn) = mean2(im_float_binned);
        proj_std(nn) = std2(im_float_binned);
    end
    fprintf('\n duration : %.1f s (%.2f min)',toc - t,(toc - t) / 60)
    t = toc;
    stat.proj_min = min(proj_min);
    stat.proj_max = max(proj_max);
    %fprintf('\npostprocessing projections:')
    toc_bytes.read_proj = tocBytes(gcp,startS);
    % Plot image statistics
    if par.visual_output
        name = 'image statistics: projections';
        his = figure('Name',name,'WindowState',window_state);

        subplot(2,2,1);
        plot(proj_min,'.')
        axis tight
        xlabel('image no.')
        title('min')

        subplot(2,2,2);
        plot(proj_max,'.')
        axis tight
        xlabel('image no.')
        title('max')

        subplot(2,2,3);
        plot(proj_mean,'.')
        axis tight
        xlabel('image no.')
        title('mean')

        subplot(2,2,4);
        plot(proj_std,'.')
        axis tight
        xlabel('image no.')
        title('std')
        drawnow
        CheckAndMakePath(fig_path)
        fig_filename = sprintf('%sfig%02u_%s.png',fig_path,his.Number,regexprep(his.Name,'\ |:','_'));
        saveas(his,fig_filename);
    end
    % Plot image statistics
    if par.visual_output && ~isempty(cur)
        name = 'image statistics: flat fields and projections';
        his = figure('Name',name,'WindowState',window_state);

        ref_ind = [cur.ref(par.ref_range).ind];
        proj_ind = [cur.proj(par.proj_range).ind];

        subplot(2,2,1);
        plot(ref_ind,flat_min,'.',proj_ind,proj_min,'.')
        legend({'ref','proj'})
        axis tight
        xlabel('image no.')
        title('min')

        subplot(2,2,2);
        plot(ref_ind,flat_max,'.',proj_ind,proj_max,'.')
        legend({'ref','proj'})
        axis tight
        xlabel('image no.')
        title('max')

        subplot(2,2,3);
        plot(ref_ind,flat_mean,'.',proj_ind,proj_mean,'.')
        legend({'ref','proj'})
        axis tight
        xlabel('image no.')
        title('mean')

        subplot(2,2,4);
        plot(ref_ind,flat_std,'.',proj_ind,proj_std,'.')
        legend({'ref','proj'})
        axis tight
        xlabel('image no.')
        title('std')

        drawnow
        CheckAndMakePath(fig_path)
        fig_filename = sprintf('%sfig%02u_%s.png',fig_path,his.Number,regexprep(his.Name,'\ |:','_'));
        saveas(his,fig_filename);
    end
    % Delete empty projections
    zz = ~projs_to_use;
    if sum(zz(:))
        fprintf('\n Deleting empty projections')
        proj(:,:,zz) = [];
    end
    if offset_shift ~= 0
        offset_shift(~projs_to_use) = [];
        x0(~projs_to_use) = [];
    end
    if vert_shift ~= 0
        vert_shift(~projs_to_use) = [];
    end
    if isfield(tomo,'vert_shift')
        if isempty(tomo.vert_shift)
            tomo.vert_shift = vert_shift;
        end
    else
        tomo.vert_shift = [];
    end
    if ~isempty(scan_position)
        scan_position(~projs_to_use) = [];
    end
    if exist('scan_position_index','var') && ~isempty(scan_position_index)
        scan_position_index(~projs_to_use) = [];
    end
    par.offset_shift_x0 = x0;
    par.offset_shift_x1 = x1;
    tomo.scan_position = scan_position;
    stat.raw_min = min(proj(:));
    stat.raw_max = max(proj(:));
    % Ring current normalization
    if par.ring_current_normalization && isempty(par.ref_path)
        proj_ind_from_filenames = proj_nums;
        proj_ind_from_log = [cur.proj(par.proj_range).ind];
        if isequal(proj_ind_from_filenames, proj_ind_from_log)
            proj_rc = [cur.proj(par.proj_range).val];
            proj_ind = [cur.proj(par.proj_range).ind];
            proj_rcm = mean(proj_rc(:));
            scale_factor = 100 ./ shiftdim(proj_rc(projs_to_use),-1);
            fprintf('\n Ring current normalization')
            %proj = fun_times(proj,scale_factor);
            proj = scale_factor .* proj;
            fprintf('\n duration : %.1f s (%.2f min)',toc - t,(toc - t) / 60)
            % Plot ring current
            if par.visual_output && exist('ref_rc','var')
                name = 'PETRA III beam current: Interpolation at image time stamps';
                if exist('hrc','var') && isvalid(hrc)
                    figure(hrc)
                else
                    hrc = figure('Name',name,'WindowState',window_state);
                end
                subplot(1,1,1);
                %ref_nums = 1:numel(ref_names(par.ref_range));
                plot(ref_ind(:),ref_rc(:),'.',proj_ind(:),proj_rc(:),'.')
                axis tight
                xlabel('image no.')
                ylabel('current / mA')
                title(name)
                legend(sprintf('flats,mean: %.2f mA',ref_rcm),sprintf('projs,mean: %.2f mA',proj_rcm))
                drawnow
                CheckAndMakePath(fig_path)
                fig_filename = sprintf('%sfig%02u_%s.png',fig_path,hrc.Number,regexprep(hrc.Name,'\ |:','_'));
                saveas(hrc,fig_filename);
            end
        else
            cprintf('Red','\n WARNING: Projections not normalized by ring current. Names read from directory and log-file are inconsistent.')
        end
    else
        cprintf('Red','\n Projections not normalized by ring current.')
    end
    stat.raw_ring_current_corrected_min = min(proj(:));
    stat.raw_ring_current_corrected_max = max(proj(:));
    num_empty = sum(~projs_to_use(:));
    par.num_proj_used = par.num_proj_used - num_empty;
    fprintf('\n crop left  min/max : %u %u',min(x0),max(x0));
    fprintf('\n crop right min/max : %u %u',min(x1),max(x1));
    fprintf('\n image shape cropbin1 : %u',im_shape_cropbin1);
    PrintVerbose(num_empty*par.delete_empty_projections,'\n discarded empty projections : %u,%.2f%%',num_empty,100*num_empty/size(proj,3))
    % if sum(num_zeros)
    %     fprintf('\n projections with zeros :')
    %     % print #zeros if not all pixels are zero
    %     for nn = 1:numel(num_zeros)
    %         if num_zeros(nn) ~= 0
    %             if isequal(num_zeros(nn),numel_im_roi_binned)
    %                 fprintf(' %u',nn)
    %             else
    %                 fprintf(' %u:%u',nn,num_zeros(nn))
    %             end
    %         end
    %     end
    % end
    fprintf('\n hot- / dark-pixel filter threshold : %f,%f',filt_pix_par.threshold_hot,filt_pix_par.threshold_dark)
    fprintf('\n global min/max of projs after filtering and binning:  %6g %6g',stat.raw_min,stat.raw_max)
    fprintf('\n global min/max of projs after dark-field correction and ring current normalization:  %6g %6g',stat.raw_ring_current_corrected_min,stat.raw_ring_current_corrected_max)
    % Save projections
    CheckAndMakePath(im_path1)
    CheckAndMakePath(im_path2)
    %% Figure: image correlation roiference
    if par.visual_output && ~sum(strcmp(image_correlation.method,{'median','mean','none'}))
        if exist('h_corr_roi','var') && isvalid(h_corr_roi)
            figure(h_corr_roi)
        else
            %h_corr_roi = figure('Name','image correlation roi','WindowState',window_state);
            figure('Name','image correlation roi','WindowState',window_state);
        end
        subplot(2,2,3)
        im = raw1;
        imsc1(im)
        hold on
        x1 = ceil(flat_corr_area1(1) / raw_bin);
        y1 = im_shape_binned2 - floor(flat_corr_area2(end) / raw_bin);
        rectangle('position',[ x1 y1 2 * flat_corr_area_width_binned 2 * flat_corr_area_height_binned],'EdgeColor','r')
        title(sprintf('proj(:,:,1)\n(vertical coordiante values flipped)'))
        colorbar
        axis equal tight

        subplot(2,2,4)
        imsc1(roi_proj(:,:,1))
        %title(sprintf('proj(:,:,1) ROI \nimage correlation area:\nunbinned (relative) coordinates\nwidth = %g %g\nheight = %g %g',image_correlation.area_width,image_correlation.area_height));
        title('proj(:,:,1) ROI')
        ylabel(sprintf('image correlation area:\nunbinned (relative) coordinates\nwidth = %g %g\nheight = %g %g',image_correlation.area_width,image_correlation.area_height));
        colorbar
        axis equal tight
        drawnow
        CheckAndMakePath(fig_path)
        fig_filename = sprintf('%sfig%02u_%s.png',fig_path,h_corr_roi.Number,regexprep(h_corr_roi.Name,'\ |:','_'));
        saveas(h_corr_roi,fig_filename);
    end
    %% Projection/flat field correlation and flat field correction %%%%%%%%
    [proj,~,toc_bytes] = proj_flat_correlation(proj,flat,image_correlation,par,write,roi_proj,roi_flat,toc_bytes);
    %[proj,corr,toc_bytes] = proj_flat_correlation(proj,flat,image_correlation,par,write,roi_proj,roi_flat,toc_bytes);
    %%%% STOP HERE TO CHECK FLATFIELD CORRELATION MAPPING %%%%%%%%%%%%%%%%%
    %%%% use 'proj_flat_sequ' to show results of the correlation
    proj_min0 = min(proj(:));
    proj_max0 = max(proj(:));
    fprintf('\n global min/max after flat-field corrected:  %6g %6g',proj_min0,proj_max0);
    write32bitTIFfromSingle(sprintf('%sproj_flatcorrected_first_%06u.tif',im_path3,1),rot90(proj(:,:,1)))
    write32bitTIFfromSingle(sprintf('%sproj_flatcorrected_last_%06u.tif',im_path3,size(proj,3)),rot90(proj(:,:,end)));
    write32bitTIFfromSingle(sprintf('%sproj_flatcorrected_mean1.tif',im_path1),squeeze(mean(proj,1)));
    write32bitTIFfromSingle(sprintf('%sproj_flatcorrected_mean2.tif',im_path2),squeeze(mean(proj,2)));
    write32bitTIFfromSingle(sprintf('%sproj_flatcorrected_mean3.tif',im_path3),rot90(squeeze(mean(proj,3))));
    %% Filter strong/full absorption (combine with iterative reco methods)
    if par.strong_abs_thresh < 1
        strong_abs_thresh = par.strong_abs_thresh;
        t = toc;
        fprintf('\n Filter flat-corrected values below %f',par.strong_abs_thresh)
        parfor nn = 1:size(proj,3)
            im = proj(:,:,nn);
            m = im < strong_abs_thresh;
            %im(m) = 0;
            if sum(m(:)) > 0
                im(m) = mean2(im(m));
                proj(:,:,nn) = im;
            end
        end
        fprintf('\n duration : %.1f s (%.2f min)',toc - t,(toc - t) / 60)
    end
    fprintf('\n sinogram size = [%g,%g,%g]',size(proj))
    %%% Angles %%
    if ~isempty(tomo.rot_angle_full_range)
        if isscalar(tomo.rot_angle_full_range)
            angles = tomo.rot_angle_full_range * (0:par.num_proj_found - 1) / par.num_proj_found;
        else
            angles = tomo.rot_angle_full_range;
        end
    else
        if isempty(imlog) || ~par.read_image_log
            if exist('cur','var') && isfield(cur,'proj') && isfield(cur.proj,'angle')
                % for KIT cam this includes missing angles
                angles = [cur.proj.angle] / 180 * pi;
                if strcmpi(cam,'kit')
                    % drop angles where projections are missing
                    angles = angles(1 + proj_nums);
                else
                    angles = angles(par.proj_range);
                end
            elseif exist(nexuslog_name{1},'file')
                %angles = s_rot.value(~boolean(stimg_key.scan.value(logpar.n_dark+1:end))) * pi / 180;
                angles = [];
                for n = 1:numel(stimg_key.scan)
                    switch numel(s_rot.value)
                        case par.num_proj_found
                            m = 1:par.num_proj_found;
                        case par.num_proj_found + par.num_ref_found
                            m =  ~boolean(stimg_key.scan(n).value(n_dark+1:end));
                        case par.num_proj_found + par.num_dark + par.num_ref_found
                            m = stimg_key.scan(n).value == 0 ;
                    end
                    angles_n = pi / 180 * s_rot.value(m);
                    angles = cat(1,angles,angles_n);
                end
                fprintf('\n angles_logged / pi: %f %f %f ... %f %f %f %f %f ',angles([1 2 3 end-4:end])/pi)
                angles = angles(par.proj_range);
            else
                num_proj = logpar.num_proj;
                switch lower(cam)
                    case 'ehd'
                        angles = tomo.rot_angle_full_range * (0:num_proj - 1) / (num_proj - 1); % EHD: ok
                    case 'kit'
                        angles = tomo.rot_angle_full_range * (0:num_proj - 1) / num_proj; % KIT: ok if logpar.projections exist
                end
            end
        end
    end
    a = sort(angles);
    da = mean(a(2:end) - a(1:end-1));
    dast = std(a(2:end) - a(1:end-1));
    fprintf('\n angle increment: mean = %g mrad = %g * pi/1000',da/1000,da/pi/1000)
    fprintf('\n angle increment: std = %g mrad = %g * pi/1000',dast/1000,dast/pi/1000)
    angles123 = a(1:3);
    angles321 = a(end-4:end);
    fprintf('\n angles / pi: %f %f %f ... %f %f %f %f %f ',angles123/pi,angles321/pi)
    fprintf('\n angles: %f * (%.1f %.1f %.1f ... %.1f %.1f %.1f %.1f %.1f)',da,angles123/da,angles321/da)
    if a(end) - a(1) > 1.1*pi
        [~,ind12] = min(abs(angles - pi));
        write32bitTIFfromSingle(sprintf('%sproj_flatcorrected_middle_%06u.tif',im_path3,ind12),fliplr(rot90(proj(:,:,ind12))))
    end
    % Figure: Angles
    if par.visual_output
        name = 'Angles';
        fa = figure('Name',name,'WindowState',window_state);
        plot(1 / pi * angles,'.')
        xlabel('projection number')
        ylabel('angle / pi')
        title(name)
        axis tight
        drawnow
        CheckAndMakePath(fig_path)
        fig_filename = sprintf('%sfig%02u_%s.png',fig_path,fa.Number,regexprep(fa.Name,'\ |:','_')) ;
        saveas(fa,fig_filename);
    end
    % drop angles where projections are empty
    angles(~projs_to_use) = [];
    %% Figure & Save: sino slice
    if par.visual_output
        if exist('h1','var') && isvalid(h1)
            figure(h1)
        else
            h1 = figure('Name','data and flat-and-dark-field correction','WindowState',window_state);
        end
        subplot(2,3,4)
        imsc1(FilterOutlier(proj(:,:,1),0.005))
        xticks([])
        yticks([])
        title(sprintf('intensity: turn 1,first proj'))
        colorbar
        axis equal tight
        subplot(2,3,5)
        if exist('scan_position_index','var') && ~isempty(scan_position_index)
            ind1 = 1:size(proj,3);
            ind1 = ind1(scan_position_index == 1);
            ind_mid = ind1(round(numel(ind1)/2));
            ind_last = ind1(end);
        else
            ind_mid = round(size(proj,3)/2);
            ind_last = size(proj,3);
        end
        imsc1(FilterOutlier(proj(:,:,ind_mid),0.005))
        xticks([])
        yticks([])
        title(sprintf('intensity: turn 1,middle proj %u',ind_mid))
        colorbar
        axis equal tight
        subplot(2,3,6)
        imsc1(FilterOutlier(proj(:,:,ind_last),0.005))
        xticks([])
        yticks([])
        title(sprintf('intensity: turn 1,last proj %u',ind_last))
        colorbar
        axis equal tight
        CheckAndMakePath(fig_path)
        fig_filename = sprintf('%sfig%02u_%s.png',fig_path,h1.Number,regexprep(h1.Name,'\ |:','_'));
        saveas(h1,fig_filename);

        h = figure('Name','Fourier transformed projection #1','WindowState',window_state);
        imf = padarray( proj(:,:,1), size(proj(:,:,1)), 'symmetric', 'post');
        imf = SubtractMean(imf);
        imf = fft2(imf);
        imf = fftshift(imf);
        imf = log(10^(-0)+abs(imf));
        imsc1(imf)
        xticks([])
        yticks([])
        title(sprintf('Fourier transformed projection'))
        colorbar
        axis equal tight
        CheckAndMakePath(fig_path)
        fig_filename = sprintf('%sfig%02u_%s.png',fig_path,h.Number,regexprep(h.Name,'\ |:','_'));
        saveas(h,fig_filename);
    end


    %% Sinogramm
    nn = round(size(proj,2) / 2);
    clear sino_mid
    if isscalar(vert_shift)
        sino_mid = squeeze(proj(:,nn,:));
    else
        y = round(nn + vert_shift);
        y(y < 1) = 1;
        y(y > size(proj,2)) = size(proj,2);
        for mm = size(proj,3):-1:1
            sino_mid(:,mm) = proj(:,y(mm),mm);
        end
    end
    sino_mean = squeeze(mean(proj,2));
    CheckAndMakePath(write.sino_path)
    filename = sprintf('%ssino_rawBin%u_middle.tif',write.sino_par_path,raw_bin);
    write32bitTIFfromSingle(filename,sino_mid);
    filename = sprintf('%ssino_rawBin%u_mean.tif',write.sino_par_path,raw_bin);
    write32bitTIFfromSingle(filename,sino_mean);
    %% Normalize sinogramm
    if par.norm_sino
        t = toc;
        fprintf('\nNormalize sinogram ')
        parfor nn = 1:size(proj,2)
            sino = proj(:,nn,:);
            sino_mean1 = mean(sino,1);
            sino_mean0 = mean(sino_mean1,3);
            m = sino_mean0 ./ sino_mean1;
            %sino_norm = bsxfun(@times,sino,m);
            sino = fun_times(sino,m);
            proj(:,nn,:) = sino;
        end
        fprintf('\n duration : %.1f s (%.2f min)',toc - t,(toc - t) / 60)
        nn = round(size(proj,2) / 2);
        if isscalar(vert_shift)
            sino_mid_norm = squeeze(proj(:,nn,:));
        else
            y = round(nn + vert_shift);
            y(y < 1) = 1;
            y(y > size(proj,2)) = size(proj,2);
            for mm = size(proj,3):-1:1
                sino_mid_norm(:,mm) = proj(:,y(mm),mm);
            end
        end
        CheckAndMakePath(write.sino_path)
        filename = sprintf('%ssino_rawBin%u_middle_normalized.tif',write.sino_par_path,raw_bin);
        write32bitTIFfromSingle(filename,sino_mid_norm);
    end
    if par.visual_output
        nn = round(size(proj,2) / 2);
        f = figure('Name','sinogram','WindowState',window_state);
        if ~par.norm_sino
            imsc1(sino_mid)
            title(sprintf('sinogram: proj(:,%u,:)',nn))
            colorbar
            axis equal tight
        else
            subplot(1,3,1)
            imsc1(sino_mid)
            title(sprintf('sinogram not normalized: proj(:,%u,:)',nn))
            colorbar
            axis equal tight
            subplot(1,3,2)
            imsc1(sino_mid_norm)
            title(sprintf('sinogram  normalized: proj(:,%u,:)',nn))
            colorbar
            axis equal tight
            subplot(1,3,3)
            imsc1(sino_mid - sino_mid_norm)
            title(sprintf('difference: map'))
            colorbar
            axis equal tight
        end
        fig_filename = sprintf('%sfig%02u_%s.png',fig_path,f.Number,regexprep(f.Name,'\ |:','_'));
        saveas(f,fig_filename);
        drawnow
        f = figure('Name','sinogram mean','WindowState',window_state);
        imsc1(sino_mean)
        title('sinogram mean')
        colorbar
        x = [-0.5,-0.25,0,0.25,0.5];
        for n = numel(x):-1:1
            xc{n} = num2str(x(n));
        end
        xt = single(normat(x) * (size(sino_mean,1) - 1) + 1);
        xticks(xt)
        xticklabels(xc)
        axis equal tight
        fig_filename = sprintf('%sfig%02u_%s.png',fig_path,f.Number,regexprep(f.Name,'\ |:','_'));
        saveas(f,fig_filename);
        drawnow
    end
    %% Ring artifact filter
    if ring_filter.apply && ring_filter.apply_before_stitching
        proj = pp_filter_ring_artefacts(ring_filter,proj,angles,par);
    end
    stat.raw_ring_current_corrected_min = min(proj(:));
    stat.raw_ring_current_corrected_max = max(proj(:));
    %% Write corrected projections
    if write.flatcor
        t = toc;
        fprintf('\nSave flat-corrected projections.')
        CheckAndMakePath(flatcor_path,write.deleteFiles,write.beamtimeID)
        parfor nn = 1:size(proj,3)
            filename = sprintf('%sproj_%s_%06u.tif',flatcor_path,scan_name,nn);
            write32bitTIFfromSingle(filename,rot90(proj(:,:,nn)));
        end
        fprintf('\n duration : %.1f (%.2f min)',toc-t,(toc-t)/60)
    end
else
    %% Read sinogram %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pp_read_sino
    %% Read flat corrected projections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pp_read_flatcor
end
%% Phase retrieval before interactive mode %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tint_phase = 0;
phase_retrieval.energy = par.energy;
phase_retrieval.sample_detector_distance = par.sample_detector_distance;
phase_retrieval.eff_pixel_size_binned = par.eff_pixel_size_binned;
if phase_retrieval.apply
    if phase_retrieval.apply_before
        %[proj,write,tomo,tint_phase] = pp_phase_retrieval(proj,phase_retrieval,tomo,write,interactive_mode);
        pp_phase_retrieval
    end
end
% %% Scaling of pixel size if MTF is wrong. Important for lateral shift and helical scans
% if ~isempty(par.pixel_scaling)
%     par.eff_pixel_size = par.pixel_scaling * par.eff_pixel_size;
%     par.eff_pixel_size_binned = par.pixel_scaling * par.eff_pixel_size_binned;
%     offset_shift = offset_shift / par.pixel_scaling;
%     scan_position = scan_position / par.pixel_scaling;
%     tomo.scan_position = scan_position;
%     fprintf('\n pixel scaling : %f ',par.pixel_scaling)
% end
%% TOMOGRAPHY: interactive mode to find rotation axis offset and tilt %%%%%
if ~par.crop_proj
    if ~isempty(tomo.rot_axis_offset_shift) && ~isscalar(tomo.rot_axis_offset_shift)
        %offset_shift = tomo.rot_axis_offset_shift / raw_bin * (0:par.num_proj_used) / par.num_proj_used;
        offset_shift = tomo.rot_axis_offset_shift / raw_bin;% * (0:par.num_proj_used) / par.num_proj_used;
        tomo.offset_shift = SubtractMean(offset_shift);
    else
        tomo.offset_shift = 1 / raw_bin * SubtractMean(offset_shift);
    end
end
[tomo.vol_shape,tomo.vol_size] = volshape_volsize(proj,tomo.vol_shape,tomo.vol_size,tomo.rot_axis_offset,verbose);
if par.virt_s_pos
    pos_s_pos_x = pos_s_pos_x_mm / 1000 / par.eff_pixel_size_binned;
    pos_s_pos_y = pos_s_pos_y_mm / 1000 / par.eff_pixel_size_binned;
    fprintf('\n pos_s_pos_x : %f mm = %f pixel',pos_s_pos_x_mm,pos_s_pos_x)
    fprintf('\n pos_s_pos_y : %f mm = %f pixel',pos_s_pos_y_mm,pos_s_pos_y)
    sh = -round(pos_s_pos_y);
    sv = -round(pos_s_pos_x);
    tomo.vol_size = tomo.vol_size + [ sh,sh,sv,sv,0,0];
    fprintf('\n tomo.vol_size shifted : ')
    fprintf(' %.1f',tomo.vol_size)
end
if ~exist('angles','var') || isempty(angles)
    d = dir([par.nexus_path filesep '*.h5']);
    nexuslog_name = [d.folder filesep d.name];
    s_rot.value = h5read(nexuslog_name,'/entry/scan/data/s_rot/value');
    par.num_dark = h5read(nexuslog_name,'/entry/scan/n_dark');
    [~,stimg_key,~,~] = pp_stimg_petra({nexuslog_name},par);
    angles = s_rot.value(~boolean(stimg_key.scan.value(par.num_dark+1:end))) * pi / 180;
end
[tomo,angles,tint,par] = interactive_mode_rot_axis(par,logpar,phase_retrieval,tomo,write,interactive_mode,proj,angles);
%% Automatic rot axis determination
tomo.angles = angles;
tomo = find_rot_axis_offset_auto(tomo,proj,par,write,interactive_mode);
%% Distortion correction
if par.distortion_correction_distance ~= 0  && ~isempty(par.distortion_correction_outer_offset)
    fprintf('\nQuadratic distortion correction')
    t = toc;
    dist_offset = par.distortion_correction_distance;
    outer_offset = par.distortion_correction_outer_offset;
    offset_diff = outer_offset - tomo.rot_axis_offset;
    exponent = par.distortion_correction_exponent;
    fprintf('\n distance: %.1f pixels',dist_offset)
    fprintf('\n rotation axis offset: %.1f pixels',tomo.rot_axis_offset)
    fprintf('\n outer rotation axis offset: %.1f pixels',outer_offset)
    fprintf('\n rotation axis offset difference: %.1f pixels',offset_diff)
    x = (1:size(proj,1)) - tomo.rot_axis_position;
    fprintf('\n orgiginal grid: x(rot axis pos) = %.1f',x(round(tomo.rot_axis_position)))
    % excentric rot axis left/right
    rao = tomo.rot_axis_offset;
    if tomo.rot_axis_offset > 0
        xq = x - 2 * offset_diff * abs(x / dist_offset).^exponent;
        fprintf('\n query grid:    xq(rot axis pos) = %.1f',xq(round(tomo.rot_axis_position)))
        fprintf('\n original grid: x(rot axis pos - dist offset) = %.1f',x(round(tomo.rot_axis_position - dist_offset)))
        fprintf('\n query grid:    xq(rot axis pos - dist offset) = %.1f',xq(round(tomo.rot_axis_position - dist_offset)))
        % Only correct from the rot axis pos to farther image edge
        xq(x>0) = x(x>0);
        % extrapolation
        o = double(ceil((x(1) - xq(1))));% max(x-xq);
    else
        xq = x + 2 * offset_diff * abs(x / dist_offset).^exponent;
        fprintf('\n query grid:    xq(rot axis pos) = %.1f',xq(round(tomo.rot_axis_position)))
        fprintf('\n orgiginal grid: x(rot axis pos + dist offset) = %.1f',x(round(tomo.rot_axis_position + dist_offset)))
        fprintf('\n query grid:    xq(rot axis pos + dist offset) = %.1f',xq(round(tomo.rot_axis_position + dist_offset)))
        % Only correct from the rot axis pos to farther image edge
        xq(x<0) = x(x<0);
        o = double(ceil((xq(end) - x(end))));% max(x-xq);
    end
    if par.visual_output
        f = figure('Name','distortion correction','WindowState',window_state);
        subplot(1,3,1)
        plot(xq - x)
        title(sprintf('displacement offset: xq - x'))
        xlabel('grid / pixel')
        ylabel('rotation axis offset difference')
        im0 = proj(:,:,1);
        drawnow
    end
    parfor nn = 1:size(proj,3)
        im = proj(:,:,nn);
        % imc = interp1(x,im,xq,'linear',1);
        if o > 0
            if rao > 0
                % extend support vector
                xo = [(x(1) + (-o:-1)),x];
                %pad image for extended vector
                imo = padarray(im,[o 0],'symmetric','pre');
            else
                %% TO BE TESTED
                xo = [x,(x(end) + (1:o))];
                imo = padarray(im,[o 0],'symmetric','post');
            end
            imc = interp1(xo,imo,xq,'linear');
        else
            imc = interp1(x,im,xq,'linear');
        end
        proj(:,:,nn) = imc;
    end
    if par.visual_output
        figure(f)
        subplot(1,3,2:3)
        imsc1(proj(:,:,1) - im0)
        axis equal tight
        xticks('auto'),yticks('auto')
        title(sprintf('difference map of proj #1: Vq - V'))
        drawnow
        fig_filename = sprintf('%sfig%02u_%s.png',fig_path,f.Number,regexprep(f.Name,'\ |:','_'));
        saveas(f,fig_filename);
    end
    fprintf('\n duration : %.1f (%.2f min)',toc-t,(toc-t)/60)
end
%% Stitch projections
if par.stitch_projections
    t = toc;
    fprintf('\nStitch projections:')
    fprintf('\n method: %s',par.stitch_method)
    num_scan_pos = max(scan_position_index);
    % LATERAL SCANNING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isscalar(offset_shift) && (num_scan_pos > 1)
        %error('Not yet implemented')
        m = angles == min(angles); % indices at same angle
        num_pos = sum(m);
        num_proj_sti = numel(angles) / num_pos;
        fprintf('\n number of scan positions incl. multi-rotation : %u',num_scan_pos)
        fprintf('\n number of lateral positions : %u',num_pos)
        fprintf('\n number projections stitched : %u',num_proj_sti)
        if round(num_proj_sti) * num_pos ~= par.num_proj_used
            error('Number of angles does not match scan positions')
        end
        angles_sti = angles(scan_position == scan_position(1)); % unique angles
        fprintf('\n number projections unstitched: %u',numel(angles))
        scan_pos_pix = scan_position - min(scan_position) + 1; % pixelwise
        fprintf('\n scan positions / pixel :')
        fprintf(' %f',scan_pos_pix(m))
        im_pos_pix = -scan_pos_pix;
        im_pos_pix = im_pos_pix - min(im_pos_pix) + 1;
        fprintf('\n image positions / pixel :')
        fprintf(' %f',im_pos_pix(m))
        pos_pix = im_pos_pix(m);
        pos_dist = pos_pix(2:end) - pos_pix(1:end-1);
        fprintf('\n distance between positions / pixel : ')
        fprintf(' %f ',pos_dist)
        if sum(pos_dist(2:end) - pos_dist(1:end-1))
            error('Unequal distance between scan positions')
        end
        % stiched image: size and preallocation
        im_shape_sti1 = floor(max(im_pos_pix) + im_shape_cropbin1-1);
        im_sti = zeros([im_shape_sti1,im_shape_binned2],'single');
        fprintf('\n image shape unstitched : %u %u',size(proj,1),size(proj,2))
        fprintf('\n image shape stitched : %u %u',size(im_sti))
        proj_sti = zeros([im_shape_sti1,im_shape_binned2,num_proj_sti],'single');
        fprintf('\n projections stitched shape: %u %u %u',size(proj_sti))
        fprintf('\n projections sitchted memory allocated : %.2f GiB',Bytes(proj_sti,3))
        ind = 1:size(proj,3);
        for nn = 1:num_proj_sti
            % indices at same angle and different positions
            m = angles == angles_sti(nn) ;
            ang_ind = ind(m);
            % stitched image
            im_sti = zeros([im_shape_sti1,im_shape_binned2],'single');

            switch lower(par.stitch_method)
                case 'step'
                    ow = par.stitch_align_overlap;
                    for pp = 1:num_pos
                        % absolute position relative to first pixel
                        x0 = pos_pix(pp);
                        x1 = pos_pix(pp) + im_shape_cropbin1 - 1;
                        x = x0:1:x1;
                        if pp < num_pos
                            x0_next = pos_pix(pp+1);
                            x1_next = pos_pix(pp+1) + im_shape_cropbin1 - 1;
                            if nn == 1
                                fprintf('\n lat pos: %u',pp)
                                fprintf('\n %8s: %6g','x0',x0)
                                fprintf('\n %8s: %6g','x1',x1)
                                fprintf('\n %8s: %6g','x0_next',x0_next)
                                fprintf('\n %8s: %6g','x1_next',x1_next)
                            end
                            if x0 < x0_next
                                % overlap = x1:x0_next;
                                % no = numel(overlap);
                                error('Not yet implemented')
                            else
                                overlap = x0:x1_next;
                                no2 = floor(numel(overlap)/2);
                                x0_o2 = overlap(no2);
                                if pp == 1
                                    x1_o2 = x1;
                                end
                                if nn == 1
                                    fprintf('\n %8s: %6g','x0_o2',x0_o2)
                                    fprintf('\n %8s: %6g','x1_o2',x1_o2)
                                end
                                xq = ceil(x0_o2):floor(x1_o2);
                                % calculate now, but used in next loop
                                x1_o2 = x0_o2 - 1;
                            end
                        else
                            x0_o2 = 1;
                            xq = ceil(x0_o2):floor(x1_o2);
                        end
                        % unstitched proj
                        v = proj(:,:,ang_ind(pp));
                        if nn == 1
                            fprintf('\n pos %u : xq = [%u %u]',pp,xq([1 end]))
                        end
                        vq = interp1(x,v,xq,'linear','extrap');
                        if x0 < x0_next
                            error('Not yet implemented')
                        else
                            r2 = ceil(0.2 * size(vq,2));
                            if pp > 1
                                m1 = mean2(vq(end-ow:end,r2:end-r2));
                                im_sti = im_sti / m2 * (m1 + m2) / 2;
                                vq = vq / m1 * (m1 + m2) / 2;
                                m2 = mean2(vq(1:ow,r2:end-r2));
                            else
                                m2 = mean2(vq(1:ow,r2:end-r2));
                            end
                        end
                        im_sti(xq,:) = vq;
                    end

                case 'linear'
                    % normalization vector
                    vec_norm = zeros([im_shape_sti1,1],'single');
                    for pp = 1:num_pos
                        % absolute position relative to first pixel
                        x0 =  pos_pix(pp);
                        x1 = pos_pix(pp) + im_shape_cropbin1 - 1;
                        xq = ceil(x0):floor(x1);
                        x = x0:1:x1;
                        % normalization vector
                        vec_norm(xq) = vec_norm(xq) + 1;
                        % unstitched proj
                        v = proj(:,:,ang_ind(pp));
                        if nn == 1
                            fprintf('\n pos %u : xq = [%u %u]',pp,xq([1 end]))
                        end
                        if mod(x0,1) == 0
                            vq = v;
                        else
                            vq = interp1(x,v,xq,'linear','extrap');
                        end
                        im_sti(xq,:) = im_sti(xq,:) + vq;
                    end
                    im_sti = im_sti ./ vec_norm;
            end
            proj_sti(:,:,nn) = im_sti;
            % Show stitched projections
            if nn == 1 || nn == num_proj_sti
                p = sprintf('%sstitched/',im_path3);
                CheckAndMakePath(p)
                write32bitTIFfromSingle(sprintf('%sproj_flatcorrected_stitched_%06u.tif',p,nn),rot90(im_sti))
                if par.visual_output
                    figure('Name','First projection stitched','WindowState',window_state);
                    switch lower(par.stitch_method)
                        case 'linear'
                            subplot(5,1,1:4)
                            imsc1(im_sti)
                            axis equal tight
                            str = ['stitched projections:' sprintf(' %u',ang_ind)];
                            title(str)
                            subplot(5,1,5)
                            plot(vec_norm)
                            axis tight
                            title('overlap positions')
                        case 'step'
                            imsc1(im_sti)
                            axis equal tight
                            str = ['stitched projections:' sprintf(' %u',ang_ind)];
                            title(str)
                    end
                end
            end
        end
        angles = angles_sti;
        proj = proj_sti;
        tomo.angles = angles;
        tomo.scan_position = [];
        clear proj_sti
    else
        % OFF-CENTERED ROTATION AXIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('\n stiching at rotation axis without interpolation.')
        fprintf('\n method: %s',par.stitch_method)
        amm = max(angles) - min(angles);
        a = sort(angles);
        da = mean(a(2:end) - a(1:end-1));
        fprintf('\n angular range before stitching: [min max diff inc] = [%f %f %f %f] degree',[min(angles),max(angles),amm]/pi*180,da*180/pi)
        if amm > 2.5 * pi
            error('\n Angular range of %f too large for stitching.',amm)
        end
        % last projection within [0,pi)
        [~,num_proj_sti] = min(abs(angles - pi));
        % number of stitched projections
        num_proj_sti = num_proj_sti - 1;
        fprintf('\n corresponding angles:\n  ')
        fprintf('%12u',1:3)
        fprintf('\n  ')
        fprintf('%12f',angles(1:3)/pi*180)
        fprintf('\n  ')
        fprintf('%12f',angles(num_proj_sti + (1:3))/pi*180)
        fprintf('\n  ')
        fprintf('%12u',num_proj_sti +(1:3))
        % rot axis pos right of center
        if tomo.rot_axis_offset >= 0
            % index range of projections to be stitched
            xl = 1:round(tomo.rot_axis_position);
            xr = 1:xl(end)-1;
            im_shape_sti1 = numel(xl) + numel(xr);
            % Preallocation
            proj_sti = zeros(im_shape_sti1,size(proj,2),num_proj_sti,'single');
            y0 = ceil(0.2 * im_shape_binned2);
            y1 = ceil(0.8 * im_shape_binned2);
            y = y0:y1;
            for nn = 1:num_proj_sti
                nn2 = mod(num_proj_sti + nn,size(proj,3)) + 1;
                im = zeros(im_shape_sti1,size(proj,2));
                % overlap region
                overlap = round(2 * tomo.rot_axis_position) - im_shape_cropbin1 : im_shape_cropbin1;
                % overlap ramp
                x = (0:1/(numel(overlap)-1):1);
                % 1D weight
                w = ones(im_shape_cropbin1,1);
                % aligment weights to adjust gray values near overlap region
                if par.stitch_align_overlap > 0
                    lo = min([floor(numel(overlap)/2),par.stitch_align_overlap]);
                    xm = round(numel(overlap)/2);
                    ml = mean2(proj(overlap(xm) + (-lo:0),y,nn));
                    mr = mean2(proj(overlap(xm) + (-lo:0),y,nn2));
                    mlr = (ml + mr) / 2;
                else
                    ml = 1;
                    mr = 1;
                    mlr = 1;
                end
                switch lower(par.stitch_method)
                    case 'step'
                        im = cat(1,mlr / ml * proj(xl,:,nn),mlr / mr * flipud(proj(xr,:,nn2)));
                        if nn == 1 || nn == round(num_proj_sti/2) || nn == num_proj_sti
                            iml = mlr / ml * proj(:,:,nn);
                            imr = flipud(mlr / mr * proj(:,:,nn2));
                        end
                    case 'linear'
                        w(overlap) = 1 - x;
                        % weighted projections
                        iml = bsxfun(@times,mlr / ml * proj(:,:,nn),w);
                        imr = flipud(bsxfun(@times,mlr / mr * proj(:,:,nn2),w));
                        % stitched projection
                        im(1:im_shape_cropbin1,:) = iml;
                        im(end - im_shape_cropbin1 + 1:end,:) = im(end - im_shape_cropbin1 + 1:end,:) + imr;
                    case 'sine'
                        w1 = w;
                        w2 = w;
                        w1(overlap) = cos(pi/2*x).^2;
                        w2(1:numel(overlap)) = sin(pi/2*x).^2;
                        % weighted projections
                        iml = bsxfun(@times,mlr / ml * proj(:,:,nn),w1);
                        imr = bsxfun(@times,mlr / mr * flipud(proj(:,:,nn2)),w2);
                        % stitched projection
                        im(1:im_shape_cropbin1,:) = iml;
                        im(end - im_shape_cropbin1 + 1:end,:) = im(end - im_shape_cropbin1 + 1:end,:) + imr;
                end
                % if tomo.filter_rotaxis
                % r = 3;
                % xr = tomo.rot_axis_position + (-2*r:2*r);
                % imroi = im(xr,:);
                %     imf = im(
                %     im(x,:) = imf;
                % end
                proj_sti(:,:,nn) = im;
                if nn == 1 || nn == round(num_proj_sti/2) || nn == num_proj_sti
                    p = sprintf('%sstitched/',im_path3);
                    CheckAndMakePath(p)
                    write32bitTIFfromSingle(sprintf('%sproj_flatcorrected_stitched_%06u.tif',p,nn),rot90(im))
                    if par.visual_output
                        %%
                        f = figure('Name',sprintf('Stitching projection %u',nn),'WindowState',window_state);
                        % proj left
                        subplot(2,2,1)
                        imsc1(iml)
                        title(sprintf('projection left %u,%f*pi rad',n,angles(nn)/pi))
                        axis equal tight
                        xticks([]),yticks([])
                        % proj right
                        subplot(2,2,2)
                        imsc1(imr);
                        title(sprintf('projection right %u,%f*pi rad',n,angles(nn)/pi))
                        axis equal tight
                        xticks([]),yticks([])
                        % proj stitched
                        subplot(2,2,3:4)
                        imsc1(im);
                        title(sprintf('projection stitched %u,%f*pi rad',n,angles(nn)/pi))
                        axis equal tight
                        xticks([]),yticks([])
                        drawnow
                        %%
                        fig_filename = sprintf('%sfig%02u_%s.png',fig_path,f.Number,regexprep(f.Name,'\ |:','_'));
                        saveas(f,fig_filename);
                    end
                end
            end
        else
            % index range of projections to be stitched
            xr = round(tomo.rot_axis_position):im_shape_cropbin1;
            xl = xr(1)+1:im_shape_cropbin1;
            im_shape_sti1 = numel(xl) + numel(xr);
            % Preallocation
            proj_sti = zeros(im_shape_sti1,size(proj,2),num_proj_sti,'single');
            for nn = 1:num_proj_sti
                if nn == 2188
                    keyboard
                end
                nn2 = mod(num_proj_sti + nn,size(proj,3)) + 1;
                im = zeros(im_shape_sti1,size(proj,2));
                % overlap region
                overlap = 1:round(2 * tomo.rot_axis_position);
                % overlap ramp
                x = (0:1/(numel(overlap)-1):1);
                switch lower(par.stitch_method)
                    case 'step'
                        im = cat(1,flipud(proj(xl,:,nn2),proj(xr,:,nn)));
                    case {'linear','sine'}
                        % 1D weight
                        w = ones(im_shape_cropbin1,1);
                        switch lower(par.stitch_method)
                            case 'linear'
                                w(overlap) =  x;
                            case 'sine'
                                w(overlap) = 0.5 - 0.5 * cos(pi*x);
                        end
                        % weighted projections
                        imr = bsxfun(@times,proj(:,:,nn),w);
                        iml = flipud(bsxfun(@times,proj(:,:,nn2),w));
                        % stitched projection
                        im(1:im_shape_cropbin1,:) = iml;
                        im(end - im_shape_cropbin1 + 1:end,:) = im(end - im_shape_cropbin1 + 1:end,:) + imr;
                end
                proj_sti(:,:,nn) = im;
            end
        end
        proj = proj_sti;
        clear proj_sti;
        %proj0 = proj;
        %proj = proj_sti;
        angles = angles(1:num_proj_sti);
        tomo.angles = angles;
        tomo.rot_axis_offset_before_stitching = tomo.rot_axis_offset;
        tomo.rot_axis_offset = 0;
        tomo.rot_axis_position_before_stitching = tomo.rot_axis_position;
        tomo.rot_axis_position = size(proj,1) / 2;
    end
    fprintf('\n duration : %.1f (%.2f min)',toc-t,(toc-t)/60)
    fprintf('\n shape of stitched projections : %u %u %u',size(proj))
    fprintf('\n memory allocated : %.2f GiB',Bytes(proj,3))
    if par.visual_output
        f = figure('Name','stitched projections','WindowState',window_state);
        % proj 1st
        subplot(3,1,1)
        imsc1(proj(:,:,1));
        title(sprintf('projection #1,%f degree',angles(1)/pi*180))
        colorbar
        axis equal tight
        xticks('auto'),yticks('auto')
        % proj middle
        n12 = round(num_proj_sti/2);
        subplot(3,1,2)
        imsc1(proj(:,:,n12));
        title(sprintf('projection #%u,%f degree',n12,angles(n12)/pi*180))
        colorbar
        axis equal tight
        xticks('auto'),yticks('auto')
        % proj last
        subplot(3,1,3)
        imsc1(proj(:,:,end));
        title(sprintf('projection #%u,%f degree',num_proj_sti,angles(end)/pi*180))
        colorbar
        axis equal tight
        xticks('auto'),yticks('auto')
        drawnow
        fig_filename = sprintf('%sfig%02u_%s.png',fig_path,f.Number,regexprep(f.Name,'\ |:','_'));
        saveas(f,fig_filename);
    end
    %% Write corrected projections
    if write.flatcor_stitched
        t = toc;
        fprintf('\nSave stitched flat-corrected projections.')
        CheckAndMakePath(flatcor_path_stitched,write.deleteFiles,write.beamtimeID)
        parfor nn = 1:size(proj,3)
            filename = sprintf('%sproj_%s_%06u.tif',flatcor_path_stitched,scan_name,nn);
            write32bitTIFfromSingle(filename,rot90(proj(:,:,nn)));
        end
        fprintf('\n duration : %.1f (%.2f min)',toc-t,(toc-t)/60)
    end
end
%% Ring artifact filter %%
if ring_filter.apply && ~ring_filter.apply_before_stitching
    proj = pp_filter_ring_artefacts(ring_filter,proj,angles,par);
    stat.raw_ring_current_corrected_min = min(proj(:));
    stat.raw_ring_current_corrected_max = max(proj(:));
end
%% Crop projections at rotation axis position %%
if par.crop_at_rot_axis
    t = toc;
    fprintf('\nCropping projections:')
    % Crop projections to avoid oversampling for scans with excentric rotation axis
    % and reconstruct WITHOUT stitching
    tomo.rot_axis_position_before_cropping = tomo.rot_axis_position;
    tomo.rot_axis_offset_before_cropping = tomo.rot_axis_offset;
    % Crop relative to rot axis position
    r = tomo.rot_axis_position / im_shape_cropbin1;
    if r < 0.5
        proj(1:floor(tomo.rot_axis_position)-1,:,:) = [];
        crop1 = im_shape_binned1 - size(proj,1);
        tomo.rot_axis_position = tomo.rot_axis_position - crop1/2;
        tomo.rot_axis_offset = tomo.rot_axis_offset - crop1/2;
    else
        proj(ceil(tomo.rot_axis_position) + 1:end,:,:) = [];
        crop1 = im_shape_binned1 - size(proj,1);
        tomo.rot_axis_position = tomo.rot_axis_position + crop1/2;
        tomo.rot_axis_offset = tomo.rot_axis_offset + crop1/2;
    end
    if isempty(tomo.vol_shape)
        tomo.vol_shape = [raw_im_shape_binned1,raw_im_shape_binned1,raw_im_shape_binned2];
    end
    fprintf('\n duration : %.1f (%.2f min)',toc-t,(toc-t)/60)
end
%% Save sinogram %%
if write.sino
    t = toc;
    fprintf('\nSave sinogram:')
    CheckAndMakePath(write.sino_path,write.deleteFiles,write.beamtimeID)
    sino_path = write.sino_path;
    parfor nn = 1:size(proj,2)
        filename = sprintf('%ssino_%s_%06u.tif',sino_path,scan_name,nn);
        sino = squeeze(proj(:,nn,:))';
        write32bitTIFfromSingle(filename,sino)
    end
    pause(0.01)
    fprintf('\n duration : %.1f s (%.2f min)',toc-t,(toc-t)/60)
end
%% Phase retrieval %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if phase_retrieval.apply
    if ~phase_retrieval.apply_before
        pp_phase_retrieval
    end
end
im_path1reco = [im_path1reco write.phase_appendix filesep];
im_path2reco = [im_path2reco write.phase_appendix filesep];
im_path3reco = [im_path3reco write.phase_appendix filesep];
%% Tomographic reco %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if tomo.run
    ttomo = toc;
    fprintf('\nTomographic reconstruction:')
    fprintf('\n method : %s',tomo.algorithm)
    fprintf('\n angle scaling : %g',tomo.angle_scaling)
    fprintf('\n angles [first last]/pi : [%g %g]',angles([1 end]) / pi)
    ro = tomo.rot_axis_offset(ceil(numel(tomo.rot_axis_offset) / 2));
    rp = tomo.rot_axis_position(ceil(numel(tomo.rot_axis_position) / 2));
    fprintf('\n rotation axis offset: %.2f',ro);
    fprintf('\n rotation axis position: %.2f',rp);
    tilt_cam = tomo.rot_axis_tilt_camera;
    tilt_lam = tomo.rot_axis_tilt_lamino;
    fprintf('\n rotation axis tilt camera: %g rad (%g deg)',tilt_cam,tilt_cam * 180 / pi)
    fprintf('\n rotation axis tilt lamino: %g rad (%g deg)',tilt_lam,tilt_lam * 180 / pi)
    if tilt_cam ~= 0
        fprintf('\n rotation axis tilt across camera: [h v] = [%.1f %.1f] pixel binned',tan(tilt_cam)*[im_shape_binned1 im_shape_binned2])
    end
    %fprintf('\n rotation axis position: %.2f',rp);
    fprintf('\n volume shape : [%g %g %g]',tomo.vol_shape)
    fprintf('\n volume size : [%g %g %g %g %g %g]',tomo.vol_size)
    vol_mem = prod(tomo.vol_shape) * 4;
    fprintf('\n volume memory : %.2f GiB',vol_mem / 1024^3)
    % Change 'reco_mode' to 'slice' if low on memory
    [mem_free,mem_avail_cpu,mem_total_cpu] = free_memory;
    fprintf('\n system memory: free,available,total : %.2f GiB,%.2f GiB,%.2f GiB',mem_free/1024^3,mem_avail_cpu/1024^3,mem_total_cpu/1024^3)
    if vol_mem > 0.9 * mem_avail_cpu
        tomo.reco_mode = 'slice';
        fprintf('\nSwitch to slice-wise reco (tomo.reco_mode = ''%s'') due to limited memory (avail : %.1f GiB,vol : %.1f GiB) .',tomo.reco_mode,mem_avail_cpu / 1024^3,vol_mem / 1024^3)
    end
    % if par.stitch_projections
    %     tomo.rot_axis_offset_before_stitching = tomo.rot_axis_offset;
    %     tomo.rot_axis_offset = 0;
    %     tomo.rot_axis_position_before_stitching = tomo.rot_axis_position;
    %     tomo.rot_axis_position = size(proj,1) / 2;
    %     rot_axis_offset_reco = 0;
    %elseif par.crop_at_rot_axis
    if par.crop_at_rot_axis
        rot_axis_offset_reco = tomo.rot_axis_position - size(proj,1) / 2;
    else
        rot_axis_offset_reco = tomo.rot_axis_offset;
    end
    % Delete redundant projection and angle
    if isequal(angles(1),angles(end))
        angles(end) = [];
        proj(:,:,end) = [];
    end
    tomo.angles = angles;

    % Filter sinogram
    if strcmpi(tomo.algorithm,'fbp')
        fprintf('\nFilter sino:')
        t2 = toc;
        filt = iradonDesignFilter(tomo.fbp_filter_type,(1 + tomo.fbp_filter_padding) * size(proj,1),tomo.fbp_filter_freq_cutoff);
        if tomo.butterworth_filter
            [b,a] = butter(tomo.butterworth_filter_order,tomo.butterworth_filter_frequ_cutoff);
            tmp = abs(freqz(b,a,round(numel(filt) / 2)));
            bw = [tmp; flipud(tmp) ];
            clear tmp
            filt = filt .* bw;
        end
        proj_shape1 = size(proj,1);
        take_neg_log = tomo.take_neg_log;
        padding = tomo.fbp_filter_padding;
        padding_method = tomo.fbp_filter_padding_method;
        switch lower(tomo.reco_mode)
            case '3d'
                fprintf('\n mode: 3d')
                startS = ticBytes(gcp);
                %mem_per_proj = 2 * 2 * 2 * 4 * prod(size(proj,1));
                % copying * complex * padding * bit detpth * (dim1 * dim2)
                %poolsize_max_fbp_filter = floor(0.9 * mem_avail_cpu / mem_per_proj);
                %fprintf('\n filter poolsize limit : %u ',poolsize_max_fbp_filter)
                %parfor (nn = 1:size(proj,3),poolsize_max_fbp_filter)
                parfor nn = 1:size(proj,3)
                    im = proj(:,:,nn);
                    im = NegLog(im,take_neg_log);
                    im = padarray(im,padding * [proj_shape1 0 0],padding_method,'post');
                    im = fft(im,[],1);
                    im = fun_times(im,filt);
                    im = ifft(im,[],1,'symmetric');
                    im = real(im);
                    im = im(1:proj_shape1,:,:);
                    proj(:,:,nn) = im;
                end
                toc_bytes.fbp_filter_sino = tocBytes(gcp,startS);
            case 'slice'
                fprintf('\n mode: slice')
                fprintf('\n slice number:\n')
                nn_count = 0;
                for nn = size(proj,2):-1:1
                    if mod(nn-1,100) == 0 || nn == size(proj,3) || nn == 1
                        nn_count = nn_count + 1;
                        fprintf(' %u',nn)
                        if ~mod(nn_count,25)
                            fprintf('\n')
                        end
                    end
                    im = proj(:,nn,:);
                    im = NegLog(im,take_neg_log);
                    im = padarray(im,padding * [proj_shape1 0 0],padding_method,'post');
                    im = fft(im,[],1);
                    im = fun_times(im,filt);
                    im = ifft(im,[],1,'symmetric');
                    im = real(im);
                    im = im(1:proj_shape1,:,:);
                    proj(:,nn,:) = im;
                end
        end
        fprintf('\n duration : %.2f min.',(toc - t2) / 60)
    else
        proj = NegLog(proj,tomo.take_neg_log);
    end
    % half weight pixel at rot axis pos as it is backprojected twice
    if par.crop_at_rot_axis
        r = tomo.rot_axis_position / im_shape_cropbin1;
        if r < 0.5
            proj(1,:,:) = 0.5 * proj(1,:,:) ;
        else
            proj(end,:,:) = 0.5 * proj(end,:,:) ;
        end
    end
    switch lower(tomo.reco_mode)
        case '3d'
            vol = zeros(tomo.vol_shape,'single');
            fprintf('\n volume memory allocated for ''3D'' mode: %.2f GiB',Bytes(vol,3))
            fprintf('\nBackproject:')
            t2 = toc;
            clearGPUs;
            vol = astra_parallel3D(tomo,permute(proj,[1 3 2]));
            fprintf('\n duration : %.2f min.',(toc - t2) / 60)
            stat.vol_min = min(vol(:));
            stat.vol_max = max(vol(:));
            %% Show orthogonal vol cuts
            if par.visual_output
                % vol z
                f = figure('Name','Volume cut z','WindowState',window_state);
                nn = round(size(vol,3) / 2);
                im = squeeze(vol(:,:,nn));
                im =  FilterOutlier(im,0.01);
                imsc(im)
                axis equal tight
                title(sprintf('vol z = %u',nn))
                colorbar
                fig_filename = sprintf('%sfig%02u_%s.png',fig_path,f.Number,regexprep(f.Name,'\ |:','_'));
                saveas(f,fig_filename);
                % vol y
                f = figure('Name','Volume cut y','WindowState',window_state);
                nn = round(size(vol,2) / 2);
                im = rot90(squeeze(vol(:,nn,:)),-2);
                im = FilterOutlier(im,0.01);
                if size(im,1) < size(im,2)
                    imsc(im)
                    title(sprintf('vol y = %u',nn))
                else
                    imsc(im')
                    title(sprintf('vol y = %u,rotated',nn))
                end
                axis equal tight
                colorbar
                fig_filename = sprintf('%sfig%02u_%s.png',fig_path,f.Number,regexprep(f.Name,'\ |:','_'));
                saveas(f,fig_filename);
                % vol x
                f = figure('Name','Volume cut x','WindowState',window_state);
                nn = round(size(vol,1) / 2);
                im = flipud(squeeze(vol(nn,:,:)));
                im =  FilterOutlier(im,0.01);
                if size(im,1) < size(im,2)
                    imsc(im)
                    title(sprintf('vol x = %u',nn))
                else
                    imsc(im')
                    title(sprintf('vol x = %u rotated',nn))
                end
                axis equal tight
                colorbar
                fig_filename = sprintf('%sfig%02u_%s.png',fig_path,f.Number,regexprep(f.Name,'\ |:','_'));
                saveas(f,fig_filename);
                drawnow
            end
            %% Save ortho slices
            fprintf('\nSave ortho slices and min/max/mean/std projections:')
            t3 = toc;
            CheckAndMakePath(write.reco_path)
            imah = @(im) (adapthisteq(normat(im)));
            % Save ortho slices x
            if ndims(vol) == 3
                nn = round(size(vol,1) / 2);
                im = squeeze(vol(nn,:,:));
                CheckAndMakePath(im_path1reco)
                filename = sprintf('%sreco_1Mid.tif',im_path1reco);
                write32bitTIFfromSingle(filename,rot90(im,-1));
                if ~anynan(im)
                    im_imah = imah(im);
                    filename = sprintf('%sreco_1MidAdaptHisteq.tif',im_path1reco);
                    write32bitTIFfromSingle(filename,rot90(im_imah,-1));
                end
                % Save ortho slices y
                CheckAndMakePath(im_path2reco)
                nn = round(size(vol,2) / 2);
                im = squeeze(vol(:,nn,:));
                filename = sprintf('%sreco_2Mid.tif',im_path2reco);
                write32bitTIFfromSingle(filename,rot90(im,-1));
                if ~anynan(im)
                    im_imah = imah(im);
                    filename = sprintf('%sreco_2MidAdaptHisteq.tif',im_path2reco);
                    write32bitTIFfromSingle(filename,rot90(im_imah,-1));
                end
            end
            % Save ortho slices z
            CheckAndMakePath(im_path3reco)
            nn = round(size(vol,3) / 2);
            im = squeeze(vol(:,:,nn));
            filename = sprintf('%sreco_3Mid.tif',im_path3reco);
            write32bitTIFfromSingle(filename,rot90(im,0));
            filename = sprintf('%sreco_3MidAdaptHisteq.tif',im_path3reco);
            write32bitTIFfromSingle(filename,rot90(imah(im),0));
            %% Save mean and min/max projections
            for h = {{'Max',@(vol,dd) max(vol,[],dd)},...
                    {'Min',@(vol,dd) min(vol,[],dd)},...
                    {'Mean',@(vol,dd) mean(vol,dd)},...
                    {'Std',@(vol,dd) std(vol,0,dd)},...
                    }
                f = h{1}{2};
                fname = h{1}{1};
                if ndims(vol) == 3
                    d0 = 1;
                else
                    d0 = 3;
                end
                for dd = d0:3
                    im = squeeze(f(vol,dd));
                    switch outputformat
                        case 'tif'
                            if phase_retrieval.apply
                                filename = sprintf('%sreco_phase%u/%s/reco_%uProj%s.tif',im_path,dd,write.phase_appendix,dd,fname);
                            else
                                filename = sprintf('%sreco%u/reco_%uProj%s.tif',im_path,dd,dd,fname);
                            end
                            write32bitTIFfromSingle(filename,rot90(im,(dd==3) - 1));
                        case {'hdf_volume','hdf_slice'}
                            filename = sprintf('%sreco%u/reco_%uProj%s.h5',im_path,dd,dd,fname);
                            if ~exist(filename,'file')
                                h5create(filename,'/volume',size(vol),'Datatype','single')
                            end
                            h5write(filename,'/volume',vol)
                    end
                end
            end
            fprintf('\n duration : %.2f min.',(toc - t3) / 60)
            %% Save volume
            if ~isfield(write,'reco_binning_factor')
                write.reco_binning_factor = 2;
            end
            if write.reco
                reco_bin = write.reco_binning_factor; % alias for readablity
                CheckAndMakePath(write.reco_path,0)
                % Save reco path to file
                filename = [userpath,filesep,'path_to_reco'];
                %filename = [getenv('HOME'),filesep,'path_to_reco'];
                fid = fopen(filename,'w');
                fprintf(fid,'%s',write.reco_path);
                fclose(fid);
                % Single precision: 32-bit float tiff
                startS = ticBytes(gcp);
                write_volume(write.float,vol,'float',write,raw_bin,phase_bin,1,0,verbose,'');
                toc_bytes.write_float = tocBytes(gcp,startS);
                % Filtered single precision: 32-bit float tiff
                write_volume(write.float_adapthisteq,vol,'float_adapthisteq',write,raw_bin,phase_bin,1,0,verbose,'');
                % Compression of dynamic range
                if write.uint8 || write.uint8_binned || write.uint16 || write.uint16_binned
                    [tlow,thigh] = compression(vol,write.compression_method,write.compression_parameter);
                else
                    tlow = 0;
                    thigh = 1;
                end
                write.tlow = tlow;
                write.thigh = thigh;
                % 16-bit tiff
                write_volume(write.uint16,vol,'uint16',write,raw_bin,phase_bin,1,0,verbose,'');
                % 8-bit tiff
                write_volume(write.uint8,vol,'uint8',write,raw_bin,phase_bin,1,0,verbose,'');
                % Bin data
                if write.float_binned || write.uint16_binned || write.uint8_binned || write.uint8_segmented
                    fprintf('\nBinning:')
                    t2 = toc;
                    vol = Binning(vol,reco_bin) / reco_bin^3;
                    fprintf('\n duration : %.2f min.',(toc - t2) / 60)
                end
                % Binned single precision: 32-bit float tiff
                write_volume(write.float_binned,vol,'float',write,raw_bin,phase_bin,reco_bin,0,verbose,'');
                % 16-bit tiff binned
                write_volume(write.uint16_binned,vol,'uint16',write,raw_bin,phase_bin,reco_bin,0,verbose,'');
                % 8-bit tiff binned
                write_volume(write.uint8_binned,vol,'uint8',write,raw_bin,phase_bin,reco_bin,0,verbose,'');
                % segmentation
                if write.uint8_segmented
                    [vol,out] = segment_volume(vol,2^10,par.visual_output,verbose);
                    write_volume(1,vol/255,'uint8',write,raw_bin,phase_bin,reco_bin,0,verbose,'_segmented');
                    save_path = write_volume(1,vol/255,'uint8',write,raw_bin,phase_bin,reco_bin,0,verbose,'_segmented');
                    save(sprintf('%ssegmentation_info.m',save_path),'out','-mat','-v7.3')
                end
            end
        case {'slice','2d'}
            %% Slicewise backprojection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Remove vertical shift for spiral CT
            if vert_shift
                t2 = toc;
                fprintf('\n Remove vertical shift:')
                parfor nn = 1:size(proj,3)
                    im = proj(:,:,nn);
                    imt = imtranslate(im,[-vert_shift(nn) 0],'linear');
                    proj(:,:,nn) = imt;
                end
                fprintf('\n duration : %.1f s (%.2f min)',toc-t2,(toc-t2)/60)
            end
            fprintf('\nBackproject and save slices:')
            %t2 = toc;
            if write.reco
                CheckAndMakePath(write.reco_path,0)
                % Save reco path to file
                filename = [userpath,filesep,'path_to_reco'];
                %filename = [getenv('HOME'),filesep,'path_to_reco'];
                fid = fopen(filename,'w');
                fprintf(fid,'%s',write.reco_path);
                fclose(fid);
            end
            % Loop over slices
            stat.vol_min = Inf;
            stat.vol_max = -Inf;
            % ASTRA parpool GPU limit
            gpu_mem_requ_per_reco = (prod(tomo.vol_shape(1:2)) + size(proj,1) * size(proj,3)) * 4;
            poolsize_max_astra = floor(par.poolsize_gpu_limit_factor * min(mem_avail_gpu) / gpu_mem_requ_per_reco);
            %% TODO: Check max poolsize
            %poolsize_max_gpu = poolsize_max_gpu * numel(par.gpu_index);
            fprintf('\n estimated GPU memory per reco : %g MiB',gpu_mem_requ_per_reco / 1024^2)
            fprintf('\n GPU poolsize limit factor : %g',par.poolsize_gpu_limit_factor)
            fprintf('\n GPU memory induced maximum poolsize : %u ',poolsize_max_astra)
            fprintf('\n Start (parallel) GPU reco: ')
            write_reco = write.reco;
            write_float =  write.float;
            %reco_path = write.reco_path;
            num_slices = size(proj,2);
            h5_filename = '';
            switch outputformat
                case 'tif'
                case 'hdf_volume'
                    h5_filename = sprintf('%s%s.h5',write.reco_path,scan_name);
            end
            %gpu_index = par.gpu_index;
            %num_gpu = numel(gpu_index);
            %poolsize_max_astra = min([poolsize_max_astra,3 *  num_gpu]);
            %parfor (nn = 1:num_slices,poolsize_max_astra)
            nn_count = 0;
            for nn = 1:num_slices
                if mod(nn-1,100) == 0 || nn == num_slices || nn == 1 || nn_count == 0
                    nn_count = nn_count + 1;
                    fprintf(' %u',nn)
                    if ~mod(nn_count,25)
                        fprintf('\n')
                    end
                end
                if ~isscalar(rot_axis_offset_reco)
                    rotation_axis_offset = rot_axis_offset_reco(nn);
                else
                    rotation_axis_offset = rot_axis_offset_reco;
                end
                % Backproject
                gpu_ind = nn;
                clearGPUs;
                vol = rot90(astra_parallel2D(tomo,permute(proj(:,nn,:),[3 1 2]),gpu_ind,rotation_axis_offset),-1);
                stat.vol_min = min(stat.vol_min,min(vol(:)));
                stat.vol_max = max(stat.vol_max,max(vol(:)));
                % Save reconstruction
                switch outputformat
                    case 'tif'
                        if write_reco
                            % Single precision: 32-bit float tiff
                            %t = toc;
                            write_volume(write_float,vol,'float',write,raw_bin,phase_bin,1,nn,0,'');
                            %fprintf('\n duration : %.2f min.',(toc - t) / 60)
                        end
                    case 'hdf_volume'
                        [s1,s2] = size(vol);
                        if ~exist(h5_filename,'file')
                            h5create(h5_filename,'/volume',[s1 s2 num_slices],'Datatype','single')
                        end
                        start = [1 1 nn];
                        count = [s1 s2 1];
                        h5write(h5_filename,'/volume',vol,start,count)
                end
            end % parfor
            fprintf(' \n tomo reco duration : %.1f s (%.2f min)',toc - ttomo,(toc - ttomo) / 60)
    end
    % Free gpu memory occupied by parpool workers
    clearGPUs
    %% Plot data transfer from/to workers
    if par.visual_output && exist('toc_bytes','var')
        bytes_sum = [0 0];
        f = figure('Name','Parallel pool data transfer during image correlation','WindowState',window_state);
        str = {};
        fn = fieldnames(toc_bytes);
        marker = 'xs*+.od';
        marker_color = 'gbrmycw';
        %y yellow m magenta	c cyan r red g green b blue	w white k black
        for nn = 1:numel(fn)
            fnn = fn{nn};
            X = toc_bytes.(fnn);
            plot(X / 1024^3,[marker(nn) marker_color(nn)])
            hold on
            bytes_sum = bytes_sum + sum(X);
            fnn = regexprep(fnn,'_',' ');
            str = cat(2,str,{sprintf('%s: to',fnn),sprintf('%s: from',fnn)});
        end
        axis tight
        title(sprintf('Data transfer of workers in parpool. Total: %.1f GiB (to),%.1f GiB (from)',bytes_sum / 1024^3))
        xlabel('worker no.')
        ylabel('transferred data / GiB')
        legend(str)
        drawnow
        pause(0.1)
        CheckAndMakePath(fig_path)
        fig_filename = sprintf('%sfig%02u_%s.png',fig_path,f.Number,regexprep(f.Name,'\ |:','_'));
        saveas(f,fig_filename);
    end
    %% Write reco log file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CheckAndMakePath(write.reco_path)
    save(sprintf('%sangles.mat',write.reco_path),'angles');
    save(sprintf('%soffset_shift.mat',write.reco_path),'offset_shift');
    if write.reco
        logfile_path = write.reco_path;
    else
        logfile_path = write.path;
    end
    CheckAndMakePath(logfile_path)
    if write.reco
        if phase_retrieval.apply
            logfile_name = sprintf('%sreco_phase_%s_rawBin%u.log',logfile_path,write.phase_appendix,raw_bin);
        else
            logfile_name = sprintf('%sreco_rawBin%u.log',logfile_path,raw_bin);
        end
        fid = fopen(logfile_name,'w');
        fprintf(fid,'scan_name : %s\n',scan_name);
        fprintf(fid,'beamtime_id : %s\n',beamtime_id);
        fprintf(fid,'scan_path : %s\n',scan_path);
        fprintf(fid,'reco_path : %s\n',write.reco_path);
        fprintf(fid,'rot_axis_offset = par.raw_bin / %u * %f\n',raw_bin,rot_axis_offset_reco);
        fprintf(fid,'rot_axis_offset_reco : %f\n',rot_axis_offset_reco);
        fprintf(fid,'eff_pixel_size : %g micron\n',par.eff_pixel_size * 1e6);
        fprintf(fid,'eff_pixel_size_binned : %g micron\n',par.eff_pixel_size_binned * 1e6);
        fprintf(fid,'energy : %g eV\n',par.energy);
        fprintf(fid,'sample_detector_distance : %f m\n',par.sample_detector_distance);
        fprintf(fid,'s_stage_z : %f \n',s_stage_z_setup);
        fprintf(fid,'s_stage_x : %f \n',s_stage_x_setup);
        fprintf(fid,'full_reconstruction_time : %.1f s\n',toc);
        fprintf(fid,'date_of_reconstruction : %s\n',datetime);
        fprintf(fid,'user : %s\n',getenv('USER'));
        fprintf(fid,'MATLAB version : %s\n',version);
        fprintf(fid,'git commit ID : %s\n',git_commit_id);
        fprintf(fid,'platform : %s\n',computer);
        fprintf(fid,'GPUs : %s',gpuDevice().Name);
        fprintf(fid,'\n\nParamters\n');
        fprintf(fid,formattedDisplayText(par));
        fprintf(fid,'\nImage correlation\n');
        fprintf(fid,formattedDisplayText(image_correlation));
        fprintf(fid,'\nImage statistics\n');
        fprintf(fid,formattedDisplayText(stat));
        % Sinogram processing
        if par.read_sino && par.filter_sino
            fprintf(fid,'\nPixel filter sino\n');
            fprintf(fid,formattedDisplayText(pixel_filter_sino));
        end
        % Phase retrieval
        if phase_retrieval.apply
            fprintf(fid,'\nPhase retrieval\n');
            fprintf(fid,formattedDisplayText(phase_retrieval));
        end
        % Ring filter
        if ring_filter.apply
            fprintf(fid,'\nRing filter\n');
            fprintf(fid,formattedDisplayText(ring_filter));
        end
        % Tomo
        fprintf(fid,'\nTomo\n');
        fprintf(fid,formattedDisplayText(tomo));
        % Write
        fprintf(fid,'\nWrite\n');
        fprintf(fid,formattedDisplayText(write));
        % Interactive
        fprintf(fid,'\nInteractive mode\n');
        fprintf(fid,formattedDisplayText(interactive_mode));
        fprintf(fid,'\ntomo.rot_axis_offset = par.raw_bin / %u * %f\n',raw_bin,rot_axis_offset_reco);
        fclose(fid);
        % End of log file
        fprintf('\n log file : \n%s',logfile_name)
        fprintf('\n reco_path : \n%s',write.reco_path)
    end
    PrintVerbose(interactive_mode.rot_axis_pos,'\nTime elapsed in interactive rotation axis centering mode: %g s (%.2f min)',tint,tint / 60);
    PrintVerbose(interactive_mode.phase_retrieval,'\nTime elapsed in interactive phase retrieval mode: %g s (%.2f min)',tint_phase,tint_phase / 60);
    fprintf('\nTime elapsed for computation: %g s (%.2f min)',toc - tint -tint_phase,(toc - tint - tint_phase) / 60);
    fprintf('\nFINISHED: %s at %s\n',scan_name,datetime)
    if isfield(par,'quick_switch') && par.quick_switch
        cprintf('Red','\nATTENTION: Quick parameter switch was turned on!\n')
    end
    % Citations
    weblink1_url = 'https://www.nature.com/articles/nprot.2014.033';
    weblink1_name = sprintf('Moosmann et al,Nat Protoc 9,294 (2014): %s',weblink1_url);
    weblink1 = sprintf('<a href = "%s">%s</a>\n',weblink1_url,weblink1_name);
    weblink2_url = 'https://doi.org/10.5281/zenodo.5118737';
    weblink2_name = sprintf('github: %s',weblink2_url);
    weblink2 = sprintf('<a href = "%s">%s</a>\n',weblink2_url,weblink2_name);
    weblink3_url = 'https://www.astra-toolbox.com';
    weblink3_name = sprintf('ASTRA toolbox: %s',weblink3_url);
    weblink3 = sprintf('<a href = "%s">%s</a>\n',weblink3_url,weblink3_name);
    fprintf('\nCitations for data reconstruction:\n');
    fprintf(weblink1);
    fprintf(weblink2);
    fprintf(weblink3);
    fprintf('\n')
    diary off
    % END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dbclear if error
end