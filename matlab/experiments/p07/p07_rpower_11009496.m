function p07_reco_loop_rpower_11009496( SUBSETS, RUN_RECO, PRINT_PARAMETERS)

% Created on 18-Jun-2021 by moosmanj

if nargin < 1
    SUBSETS = [];
end
if nargin < 2
    RUN_RECO = 0;
end
if nargin < 3
    PRINT_PARAMETERS = '';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFAULT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% !!! QUICK SWITCH TO ALTERNATIVE SET OF PARAMETERS !!!
% !!! OVERWRITES PARAMETERS BELOW QUICK SWITCH SECTION !!!
% Just copy parameter and set quick switch to 1
par.quick_switch = 0;
% END OF QUICK SWITCH TO ALTERNATIVE SET OF PARAMETERS %%%%%%%%%%%%%%%%%%%%

pp_parameter_switch % DO NOT DELETE THIS LINE

%%% SCAN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%l%%%%%%%%%%%
par.scan_path = pwd; % string/pwd. pwd: change to directory of the scan to be reconstructed, string: absolute scan path
par.ref_path = {}; % cell of strings. Additonal data sets to be included for the correlation of projections and reference images
par.read_flatcor = 0; % read preprocessed flatfield-corrected projections. CHECK if negative log has to be taken!
par.read_flatcor_path = ''; % subfolder of 'flat_corrected' containing projections
par.read_flatcor_trafo = @(im) im; %fliplr( im ); % anonymous function applied to the image which is read e.g. @(x) rot90(x)
par.read_sino = 0; % read preprocessed sinograms. CHECK if negative log has to be taken!
par.read_sino_folder = ''; % subfolder to scan path
par.read_sino_trafo = @(x) (x);%rot90(x); % anonymous function applied to the image which is read e.g. @(x) rot90(x)
par.energy = []; % eV! if empty: read from log file (log file values can be ambiguous or even missing sometimes)
par.sample_detector_distance = []; % in m. if empty: read from log file
par.eff_pixel_size = []; %1.07e-6; % in m. if empty: read from log lfile. effective pixel size =  detector pixel size / magnification
par.pixel_scaling = []; % to account for mismatch of eff_pixel_size with, ONLY APPLIED BEFORE TOMOGRAPHIC RECONSTRUCTION, HAS TO BE CHANGED!
par.read_image_log = 0; % bool, default: 0. Read metadata from image log instead hdf5, if image log exists
%%% PREPROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.raw_bin = 3; % projection binning factor: integer
par.raw_roi = []; % vertical and/or horizontal ROI; (1,1) coordinate = top left pixel; supports absolute, relative, negative, and mixed indexing.
% []: use full image;
% [y0 y1]: vertical ROI, skips first raw_roi(1)-1 lines, reads until raw_roi(2); if raw_roi(2) < 0 reads until end - |raw_roi(2)|; relative indexing similar.
% [y0 y1 x0 x1]: vertical + horzontal ROI, each ROI as above
% if scalar = 0: auto roi, selects vertical ROI automatically. Use only for DCM. Not working for *.raw data where images are flipped and DMM data.
% if scalar < 1: Threshold is set as min(proj(:,:,[1 end])) + abs(raw_roi)*median(dark(:)). raw_roi=-1 defaults to min(proj(:,:,[1 end])) + 4*median(dark(:))
par.im_trafo = '';% 'rot90(im,1)'; % string to be evaluated after reading data in the case the image is flipped/rotated/etc due to changes at the beamline, e.g. 'rot90(im)'
% STITCHING/CROPPING only for scans without lateral movment. Legacy support
par.crop_at_rot_axis = 0; % for recos of scans with excentric rotation axis but WITHOUT projection stitching
par.stitch_projections = 0; % for 2 pi cans: stitch projection at rotation axis position. Recommended with phase retrieval to reduce artefacts. Standard absorption contrast data should work well without stitching. Subpixel stitching not supported (non-integer rotation axis position is rounded, less/no binning before reconstruction can be used to improve precision).
par.stitch_method = 'sine'; 'step';'linear'; %  ! CHECK correlation area !
% 'step' : no interpolation, use step function
% 'linear' : linear interpolation of overlap region
% 'sine' : sinusoidal interpolation of overlap region
par.proj_range = []; % range of projections to be used (from all found, if empty or 1: all, if scalar: stride, if range: start:incr:end
par.ref_range = []; % range of flat fields to be used (from all found), if empty or 1: all. if scalar: stride, if range: start:incr:end
par.wo_crop = 0; % Do not crop images to account for random lateral shift
par.virt_s_pos = 0; % Correct sample position in reconsructed volume if virtual sample position motors are used
pixel_filter_threshold_dark = [0.01 0.005]; % Dark fields: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_threshold_flat = [0.02 0.005]; % Flat fields: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_threshold_proj = [0.02 0.005]; % Raw projection: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_radius = [5 5]; % Increase only if blobs of zeros or other artefacts are expected. Can increase processing time heavily.
par.ring_current_normalization = 1; % normalize flat fields and projections by ring current
image_correlation.method = 'ssim-ml';'entropy';'median';'none';'ssim';'ssim-g';'std';'cov';'corr';'diff1-l1';'diff1-l2';'diff2-l1';'diff2-l2';'cross-entropy-12';'cross-entropy-21';'cross-entropy-x';
% Correlation of projections and flat fields. Essential for DCM data. Typically improves reconstruction quality of DMM data, too.
% Available methods ('ssim-ml'/'entropy' usually work best):
% 'none' : no correlation, for DPC
% 'median' or '': use median flat
% 'ssim-ml' : Matlab's structural similarity index (SSIM), includes Gaussian smoothing
% 'entropy' : entropy measure of proj over flat, usually similar result to SSIM, but faster
% 'ssim' : own implementation of SSIM without smoothing, usually worse, but sometimes better than 'ssim-ml'
% 'ssim-g' : 'ssim' with smoothing (Gaussian blurring)
% 'cov' : cross covariance
% 'corr' : cross correlation = normalized cross covariance
% 'std' : standard deviation of proj over flat
% 'diff1/2-l1/2': L1/L2-norm of anisotropic (diff1-l*) or isotropic (diff2-l*) difference of projections and flat fields
% 'cross-entropy-*' : asymmetric (12,21) and symmetric (x) cross entropy
image_correlation.force_calc = 0; % bool. force compuation of correlation even though a (previously computed) corrlation matrix exists
image_correlation.num_flats = 9; % number of best maching flat fields used for correction
image_correlation.area_width = [1 100];%[-100 1];% correlation area: index vector or relative/absolute position of [first pix, last pix], negative indexing is supported
image_correlation.area_height = [0.2 0.8]; % correlation area: index vector or relative/absolute position of [first pix, last pix]
image_correlation.filter = 1; % bool, filter ROI before correlation
image_correlation.filter_type = 'median'; % string. correlation ROI filter type, currently only 'median' is implemnted
image_correlation.filter_parameter = {[3 3], 'symmetric'}; % cell. filter paramaters to be parsed with {:}
% 'median' : using medfilt2, parameters: {[M N]-neighboorhood, 'padding'}
% 'wiener' : using wiender2, parameters: {[M N]-neighboorhood}
ring_filter.apply = 0; % ring artifact filter (use only for scans without lateral sample movement)
ring_filter.apply_before_stitching = 1; % ! Consider when phase retrieval is applied !
ring_filter.method = 'jm'; 'wavelet-fft';
ring_filter.waveletfft_dec_levels = 1:6; % decomposition levels for 'wavelet-fft'
ring_filter.waveletfft_wname = 'db7';'db25';'db30'; % wavelet type, see 'FilterStripesCombinedWaveletFFT' or 'waveinfo'
ring_filter.waveletfft_sigma = 3; %  suppression factor for 'wavelet-fft'
ring_filter.jm_median_width = 11; % multiple widths are applied consecutively, eg [3 11 21 31 39];
par.strong_abs_thresh = 1; % if 1: does nothing, if < 1: flat-corrected values below threshold are set to one. Try with algebratic reco techniques.
par.norm_sino = 0; % not recommended, can introduce severe artifacts, but sometimes improves quality
%%% PHASE RETRIEVAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phase_retrieval.apply = 0; % See 'PhaseFilter' for detailed description of parameters !
phase_retrieval.apply_before = 0; % before stitching, interactive mode, etc. For phase-contrast data with an excentric rotation axis phase retrieval should be done afterwards. To find the rotataion axis position use this option in a first run, and then turn it of afterwards.
phase_retrieval.post_binning_factor = 1; % Binning factor after phase retrieval, but before tomographic reconstruction
phase_retrieval.method = 'tie';'tieNLO_Schwinger';'dpc';'tie';'qp';'qpcut'; %'qp' 'ctf' 'tie' 'qp2' 'qpcut'
% Interactive phase retrieval not supported for method 'tieNLO_Schwinger'
phase_retrieval.reg_par = 1.1; % regularization parameter. larger values tend to blurrier images. smaller values tend to original data.
phase_retrieval.bin_filt = 0.1; % threshold for quasiparticle retrieval 'qp', 'qp2'
phase_retrieval.cutoff_frequ = 2 * pi; % in radian. frequency cutoff in Fourier space for 'qpcut' phase retrieval
phase_retrieval.padding = 1; % padding of intensities before phase retrieval, 0: no padding
phase_retrieval.tieNLO_Schwinger.sn = 10; % Schwinger regularization: points of support
phase_retrieval.tieNLO_Schwinger.smax = 10; % Schwinger regularization: maximumg support range
phase_retrieval.dpc_steps = 5;
phase_retrieval.dpc_bin = 4;
%%% TOMOGRAPHY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tomo.run = 1; % run tomographic reconstruction
tomo.run_interactive_mode = 1; % if tomo.run = 0, use to determine rot axis positions without processing the full tomogram;
tomo.reco_mode = '3D';'slice';  % slice-wise or full 3D backprojection. 'slice': volume must be centered at origin & no support of rotation axis tilt, reco binning, save compressed
tomo.vol_size = [];%[-1.5 1.5 -1.5 1.5 -0.5 0.5];% 6-component vector [xmin xmax ymin ymax zmin zmax], for excentric rot axis pos / extended FoV;. if empty, volume is centerd within tomo.vol_shape. unit voxel size is assumed. if smaller than 10 values are interpreted as relative size w.r.t. the detector size. Take care bout minus signs! Note that if empty vol_size is dependent on the rotation axis position.
tomo.vol_shape = []; %[1 1 1] shape (# voxels) of reconstruction volume. used for excentric rot axis pos. if empty, inferred from 'tomo.vol_size'. in absolute numbers of voxels or in relative number w.r.t. the default volume which is given by the detector width and height.
tomo.rot_angle_full_range = [];% 2 * pi ;[]; % in radians. if []: full angle of rotation including additional increment, or array of angles. if empty full rotation angles is determined automatically to pi or 2 pi
tomo.rot_angle_offset = pi; % global rotation of reconstructed volume
tomo.rot_axis_offset = []; % rotation axis offset w.r.t to the image center. Assuming the rotation axis position to be centered in the FOV for standard scan, the offset should be close to zero.
tomo.rot_axis_position = []; % if empty use automatic computation. EITHER OFFSET OR posITION MUST BE EMPTY. YOU MUST NOT USE BOTH!
tomo.rot_axis_offset_shift = []; %[]; % absolute lateral movement in pixels during fly-shift-scan, overwrite lateral shift read out from hdf5 log
tomo.flip_scan_position = 0; % for debugging
tomo.rot_axis_tilt_camera = 0; % in rad. camera tilt w.r.t rotation axis.
tomo.rot_axis_tilt_lamino = 0; % in rad. lamino tilt w.r.t beam.
tomo.rot_axis_corr_area1 = []; % ROI to correlate projections at angles 0 & pi. Use [0.75 1] or so for scans with an excentric rotation axis
tomo.rot_axis_corr_area2 = []; % ROI to correlate projections at angles 0 & pi
tomo.fbp_filter_type = 'Ram-Lak';'linear'; % see iradonDesignFilter for more options. Ram-Lak according to Kak/Slaney
tomo.fbp_filter_freq_cutoff = 1; % Cut-off frequency in Fourier space of the above FBP filter
tomo.fbp_filter_padding = 1; % symmetric padding for consistent boundary conditions, 0: no padding
tomo.fbp_filter_padding_method = 'symmetric';
tomo.butterworth_filter = 0; % use butterworth filter in addition to FBP filter
tomo.butterworth_filter_order = 1;
tomo.butterworth_filter_frequ_cutoff = 0.9;
tomo.astra_pixel_size = 1; % detector pixel size for reconstruction: if different from one 'tomo.vol_size' must to be ajusted, too!
tomo.take_neg_log = []; % take negative logarithm. if empty, use 1 for attenuation contrast, 0 for phase contrast
tomo.algorithm =  'fbp';'cgls';'sirt';'sart';'em';'fbp-astra'; % SART/EM only work for 3D reco mode
tomo.iterations = 50; % for iterateive algorithms: 'sirt', 'cgls', 'sart', 'em'
tomo.MinConstraint = []; % sirt3D/sirt2d/sart2d only. If specified, all values below MinConstraint will be set to MinConstraint. This can be used to enforce non-negative reconstructions, for example.
tomo.MaxConstraint = []; % sirt3D/sirt2d/sart2d only. If specified, all values above MaxConstraint will be set to MaxConstraint.
%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write.path = ''; %'/gpfs/petra3/scratch/moosmanj'; % absolute path were output data will be stored. !!overwrites the write.to_scratch flag. if empty uses the beamtime directory and either 'processed' or 'scratch_cc'
write.to_scratch = 0; % write to 'scratch_cc' instead of 'processed'
write.deleteFiles = 0; % delete files already existing in output folders. Useful if number or names of files differ when reprocessing.
write.beamtimeID = ''; % string (regexp),typically beamtime ID, mandatory if 'write.deleteFiles' is true (safety check)
write.scan_name_appendix = ''; % appendix to the output folder name which defaults to the scan name
write.parfolder = '';% sprintf( '%s%04u, tomo.algorithm, tomo.iterations); '';% parent folder to 'reco', 'sino', 'phase', and 'flat_corrected'
write.subfolder_flatcor = ''; % subfolder in 'flat_corrected'
write.subfolder_phase_map = ''; % subfolder in 'phase_map'
write.subfolder_sino = ''; % subfolder in 'sino'
write.subfolder_reco = ''; % subfolder in 'reco'
write.flatcor = 0; % save preprocessed flat corrected projections
write.phase_map = 0; % save phase maps (if phase retrieval is not 0)
write.sino = 0; % save sinograms (after preprocessing & before FBP filtering and phase retrieval)
write.phase_sino = 0; % save sinograms of phase maps
write.reco = 1; % save reconstructed slices (if tomo.run=1)
write.float = 1; % single precision (32-bit float) tiff
write.uint16 = 0; % save 16bit unsigned integer tiff using 'write.compression_method'
write.uint8 = 0; % save binned 8bit unsigned integer tiff using 'write.compression_method'
% Optionally save binned reconstructions, only works in '3D' reco_mode
write.float_binned = 0; % save binned single precision (32-bit float) tiff
write.uint16_binned = 0; % save binned 16bit unsigned integer tiff using 'write.compression_method'
write.uint8_binned = 0; % save binned 8bit unsigned integer tiff using 'wwrite.compression_method'
write.reco_binning_factor = 2; % IF BINNED VOLUMES ARE SAVED: binning factor of reconstructed volume
write.compression_method = 'outlier';'histo';'full'; 'std'; 'threshold'; % method to compression dynamic range into [0, 1]
write.compression_parameter = [0.02 0.02]; % compression-method specific parameters
% dynamic range is compressed s.t. new dynamic range assumes
% 'outlier' : [LOW, HIGH] = write.compression_parameter, eg. [0.01 0.03], outlier given in percent, if scalear LOW = HIGH.
% 'full' : full dynamic range is used
% 'threshold' : [LOW HIGH] = write.compression_parameter, eg. [-0.01 1]
% 'std' : NUM = write.compression_parameter, mean +/- NUM*std, dynamic range is rescaled to within -/+ NUM standard deviations around the mean value
% 'histo' : [LOW HIGH] = write.compression_parameter (100*LOW)% and (100*HIGH)% of the original histogram, e.g. [0.02 0.02]
write.uint8_segmented = 0; % experimental: threshold segmentaion for histograms with 2 distinct peaks: __/\_/\__
%%% INTERACTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.visual_output = 1; % show images and plots during reconstruction
interactive_mode.rot_axis_pos = 0; % reconstruct slices with dif+ferent rotation axis offsets
interactive_mode.rot_axis_pos_default_search_range = []; % if empty: asks for search range when entering interactive mode
interactive_mode.rot_axis_tilt = 0; % reconstruct slices with different offset AND tilts of the rotation axis
interactive_mode.rot_axis_tilt_default_search_range = []; % if empty: asks for search range when entering interactive mode
interactive_mode.lamino = 1; % find laminography tilt instead camera tilt
interactive_mode.angles = 0; % reconstruct slices with different scalings of angles
interactive_mode.angle_scaling_default_search_range = []; % if empty: use a variaton of -/+5 * (angle increment / maximum angle)
interactive_mode.slice_number = 0.5; % default slice number. if in [0,1): relative, if in (1, N]: absolute
interactive_mode.phase_retrieval = 1; % Interactive retrieval to determine regularization parameter
interactive_mode.phase_retrieval_default_search_range = []; % if empty: asks for search range when entering interactive mode, otherwise directly start with given search range
%%% HARDWARE / SOFTWARE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.use_cluster = 0; % if available: on MAXWELL nodes disp/nova/wga/wgs cluster computation can be used. Recommended only for large data sets since parpool creation and data transfer implies a lot of overhead.
par.poolsize = 0.9; % number of workers used in a local parallel pool. if 0: use current config. if >= 1: absolute number. if 0 < poolsize < 1: relative amount of all cores to be used. if SLURM scheduling is available, a default number of workers is used.
par.poolsize_gpu_limit_factor = 0.7; % Relative amount of GPU memory used for preprocessing during parloop. High values speed up Proprocessing, but increases out-of-memory failure
tomo.astra_link_data = 1; % ASTRA data objects become references to Matlab arrays. Reduces memory issues.
tomo.astra_gpu_index = []; % GPU Device index to use, Matlab notation: index starts from 1. default: [], uses all
par.gpu_index = tomo.astra_gpu_index;
par.use_gpu_in_parfor = 1;
phase_retrieval.use_parpool = 0; % bool. Disable parpool when out-of-memory error occurs during phase retrieval.
par.window_state = 'minimized';'normal';'maximized'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set default. Defines default parameter set to be used.
SET_DEFAULT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETER / DATA SETS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_path = '/asap3/petra3/gpfs/p07/2021/data/11009496/raw/';

% 80 mm
tomo.rot_axis_offset = -1.7 * 3 / par.raw_bin;
par.scan_path = [raw_path 'rpower001_fk01']; ADD

% 630 mm
tomo.rot_axis_offset = 1.60 * 3 / par.raw_bin;
par.scan_path = [raw_path 'rpower002_fk01_630mm']; ADD

par.raw_bin = 2;
tomo.rot_axis_offset = 2.9 * 3 / par.raw_bin;
par.scan_path = [raw_path 'rpower003_fk01_300mm']; ADD

par.raw_bin = 3;
tomo.rot_axis_offset = 5.1 * 3 / par.raw_bin;
par.scan_path = [raw_path 'rpower004_cr94']; ADD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.ring_current_normalization = 0;
% PETRA stability problems, beam current info lost after half of the scan
par.scan_path = [raw_path 'rpower005_cmol05'];
ADD
par.ring_current_normalization = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.scan_path = [raw_path 'rpower006_cmol145']; 

par.raw_bin = 2; 
tomo.rot_axis_offset = 5.1 * 3 / par.raw_bin;
ADD

par.raw_bin = 3;
tomo.rot_axis_offset = 5.1 * 3 / par.raw_bin;
ADD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

par.scan_path = [raw_path 'rpower006_cmol145_b']; 

par.raw_bin = 2; 
tomo.rot_axis_offset = 5.1 * 3 / par.raw_bin;
ADD

par.raw_bin = 3;
tomo.rot_axis_offset = 5.1 * 3 / par.raw_bin;
ADD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pixel scaling determined from 0.5 aperture scan
par.pixel_scaling = 0.9982;
par.scan_path = [raw_path 'rpower007_cmol145_a_oaperture0p5']; 

par.raw_bin = 2; 
tomo.rot_axis_offset = 5.1 * 3 / par.raw_bin;
ADD

par.scan_path = [raw_path 'rpower008_cmol145_a_oaperture0p5_singleFOV_proj20k']; 

par.raw_bin = 2; 
tomo.rot_axis_offset = 5.1 * 3 / par.raw_bin;
ADD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.pixel_scaling = 1;
interactive_mode.rot_axis_pos = 0;
par.raw_roi = [0.2 0.8];
par.proj_range = 2;
par.raw_bin = 4; 
tomo.rot_axis_offset = 17.0 * 3 / par.raw_bin;

par.scan_path = [raw_path 'rpower009_cm03_a']; ADD
par.scan_path = [raw_path 'rpower009_cm03_b']; ADD
par.scan_path = [raw_path 'rpower009_cm03_c']; ADD
par.scan_path = [raw_path 'rpower009_cm03_d']; ADD
par.scan_path = [raw_path 'rpower009_cm03_e']; ADD
par.scan_path = [raw_path 'rpower009_cm03_f']; ADD

%% CHECK axis
tomo.rot_axis_offset = (17.0 * 3 +  500 / 1.127298)/ par.raw_bin;
par.scan_path = [raw_path 'rpower010_bb02_a']; ADD
par.scan_path = [raw_path 'rpower010_bb02_b']; ADD

tomo.rot_axis_offset = 17.0 * 3 / par.raw_bin;
par.scan_path = [raw_path 'rpower012_oc125_a']; ADD
par.scan_path = [raw_path 'rpower012_oc125_b']; ADD
par.scan_path = [raw_path 'rpower013_en006_a']; ADD
par.scan_path = [raw_path 'rpower013_en006_b']; ADD
par.scan_path = [raw_path 'rpower013_en006_c']; ADD
par.scan_path = [raw_path 'rpower014_ra2011p']; ADD
par.scan_path = [raw_path 'rpower015_ra2011m']; ADD
par.scan_path = [raw_path 'rpower016_ra2026_a']; ADD
par.scan_path = [raw_path 'rpower016_ra2026_b']; ADD
% 2 x 1
par.scan_path = [raw_path 'rpower017_cmol3']; ADD
% 2 x 2
par.scan_path = [raw_path 'rpower018_ra2037_a']; ADD
par.scan_path = [raw_path 'rpower018_ra2037_b']; ADD
% 2 x 1
par.scan_path = [raw_path 'rpower019_tai006']; ADD
% 3 x 2
par.scan_path = [raw_path 'rpower020_cmol17_a']; ADD
par.scan_path = [raw_path 'rpower020_cmol17_b']; ADD
% 3 x 1
par.scan_path = [raw_path 'rpower021_tai007']; ADD
% 2 x 1, wrong name, real name: tai008
par.scan_path = [raw_path 'rpower022_ra2004']; ADD
% 2 x 1
par.scan_path = [raw_path 'rpower026_ra2029b']; ADD
% 2 x 1
par.scan_path = [raw_path 'rpower027_ra2004']; ADD
% 1 x 2
par.scan_path = [raw_path 'rpower028_csv041_a']; ADD
par.scan_path = [raw_path 'rpower028_csv041_b']; ADD

interactive_mode.rot_axis_pos = 0;
phase_retrieval.apply = 0;
par.wo_crop = 1;
par.pixel_scaling = 1;
par.raw_roi = [];
par.proj_range = 1;
image_correlation.method = 'ssim-ml';
image_correlation.force_calc = 0;
image_correlation.num_flats = 12; 
image_correlation.area_width = [1 100];
image_correlation.area_height = [0.6 0.9];

% 1 x 1
par.scan_path = [raw_path 'rpower025_cm02a']; 
par.raw_bin = 2; tomo.rot_axis_offset = 20.275 * 2 / par.raw_bin; ADD;
par.raw_bin = 3; tomo.rot_axis_offset = 20.275 * 2 / par.raw_bin; ADD;

% 1 x 1
par.scan_path = [raw_path 'rpower029_oc141'];
par.raw_bin = 2; tomo.rot_axis_offset = 34.5 * 2 / par.raw_bin; ADD;
par.raw_bin = 3; tomo.rot_axis_offset = 34.5 * 2 / par.raw_bin; ADD;

% 1 x 1
par.scan_path = [raw_path 'rpower030_fk02'];
par.raw_bin = 2; tomo.rot_axis_offset = 35.3 * 2 / par.raw_bin; ADD;
par.raw_bin = 3; tomo.rot_axis_offset = 35.3 * 2 / par.raw_bin; ADD;

% 1 x 1
par.scan_path = [raw_path 'rpower031_bb010'];
par.raw_bin = 2; tomo.rot_axis_offset = 36.8 * 2 / par.raw_bin; ADD;
par.raw_bin = 3; tomo.rot_axis_offset = 36.8 * 2 / par.raw_bin; ADD;

% 1 x 1
par.scan_path = [raw_path 'rpower032_cm01']; 
par.raw_bin = 2; tomo.rot_axis_offset = 37.6 * 2 / par.raw_bin; ADD;
par.raw_bin = 3; tomo.rot_axis_offset = 37.6 * 2 / par.raw_bin; ADD;

% 1 x 1
par.scan_path = [raw_path 'rpower033_bb12']; 
par.raw_bin = 2; tomo.rot_axis_offset = 38.8 * 2 / par.raw_bin; ADD;
par.raw_bin = 3; tomo.rot_axis_offset = 38.8 * 2 / par.raw_bin; ADD;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phase_retrieval.apply = 1;
interactive_mode.phase_retrieval = 0;
phase_retrieval.use_parpool = 1;
phase_retrieval.reg_par = 2.2;

% 1 x 1
par.scan_path = [raw_path 'rpower025_cm02a']; 
par.raw_bin = 2; tomo.rot_axis_offset = 20.275 * 2 / par.raw_bin; ADD;
par.raw_bin = 3; tomo.rot_axis_offset = 20.275 * 2 / par.raw_bin; ADD;

% 1 x 1
par.scan_path = [raw_path 'rpower029_oc141'];
par.raw_bin = 2; tomo.rot_axis_offset = 34.5 * 2 / par.raw_bin; ADD;
par.raw_bin = 3; tomo.rot_axis_offset = 34.5 * 2 / par.raw_bin; ADD;

% 1 x 1
par.scan_path = [raw_path 'rpower030_fk02'];
par.raw_bin = 2; tomo.rot_axis_offset = 35.3 * 2 / par.raw_bin; ADD;
par.raw_bin = 3; tomo.rot_axis_offset = 35.3 * 2 / par.raw_bin; ADD;

% 1 x 1
par.scan_path = [raw_path 'rpower031_bb010'];
par.raw_bin = 2; tomo.rot_axis_offset = 36.8 * 2 / par.raw_bin; ADD;
par.raw_bin = 3; tomo.rot_axis_offset = 36.8 * 2 / par.raw_bin; ADD;

% 1 x 1
par.scan_path = [raw_path 'rpower032_cm01']; 
par.raw_bin = 2; tomo.rot_axis_offset = 37.6 * 2 / par.raw_bin; ADD;
par.raw_bin = 3; tomo.rot_axis_offset = 37.6 * 2 / par.raw_bin; ADD;

% 1 x 1
par.scan_path = [raw_path 'rpower033_bb12']; 
par.raw_bin = 2; tomo.rot_axis_offset = 38.8 * 2 / par.raw_bin; ADD;
par.raw_bin = 3; tomo.rot_axis_offset = 38.8 * 2 / par.raw_bin; ADD;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.raw_bin = 2;
par.raw_roi = [0.3 0.6];
interactive_mode.phase_retrieval = 0;
phase_retrieval.apply = 1;
phase_retrieval.reg_par = 2.2;
par.scan_path = [raw_path 'rpower029_oc141'];
tomo.rot_axis_offset = 34.5 * 2 / par.raw_bin;
phase_retrieval.reg_par = 2.2; % regularization parameter. larger values tend to blurrier images. smaller values tend to original data.
phase_retrieval.bin_filt = 0.1; % threshold for quasiparticle retrieval 'qp', 'qp2'
phase_retrieval.cutoff_frequ = 2 * pi; % in radian. frequency cutoff in Fourier space for 'qpcut' phase retrieval
phase_retrieval.padding = 1; % padding of intensities before phase retrieval, 0: no padding
phase_retrieval.tieNLO_Schwinger.sn = 10; % Schwinger regularization: points of support
phase_retrieval.tieNLO_Schwinger.smax = 10; % Schwinger regularization: maximumg support range

phase_retrieval.method = 'tie';ADD
phase_retrieval.method = 'qp'; ADD
phase_retrieval.method = 'tieNLO_Schwinger'; ADD


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %par.scan_path = [raw_path 'rpower007_cmol145_a_oaperture0p5']; 
% par.scan_path = [raw_path 'rpower009_cm03_c'];
% par.raw_roi = [0.1 0.9];
% par.raw_bin = 4;
% par.ref_range = 4;
% par.wo_crop = 0;
% tomo.rot_axis_offset = 5.1 * 3 / par.raw_bin;
% image_correlation.method = 'median';
% 
% 
% par.pixel_scaling = 0.9972;
% write.subfolder_reco = sprintf('pixscal%.4f',par.pixel_scaling);
% ADD
% 
% par.pixel_scaling = 0.9974;
% write.subfolder_reco = sprintf('pixscal%.4f',par.pixel_scaling);
% ADD
% 
% par.pixel_scaling = 0.9976;
% write.subfolder_reco = sprintf('pixscal%.4f',par.pixel_scaling);
% ADD
% 
% par.pixel_scaling = 0.9978;
% write.subfolder_reco = sprintf('pixscal%.4f',par.pixel_scaling);
% ADD
% 
% par.pixel_scaling = 0.9980;
% write.subfolder_reco = sprintf('pixscal%.4f',par.pixel_scaling);
% ADD
% 
% par.pixel_scaling = 0.9982;
% write.subfolder_reco = sprintf('pixscal%.4f',par.pixel_scaling);
% ADD
% 
% par.pixel_scaling = 0.9984;
% write.subfolder_reco = sprintf('pixscal%.4f',par.pixel_scaling);
% ADD
% 
% par.pixel_scaling = 0.9986;
% write.subfolder_reco = sprintf('pixscal%.4f',par.pixel_scaling);
% ADD
% 
% par.pixel_scaling = 0.9988;
% write.subfolder_reco = sprintf('pixscal%.4f',par.pixel_scaling);
% ADD
% 
% par.pixel_scaling = 0.9990;
% write.subfolder_reco = sprintf('pixscal%.4f',par.pixel_scaling);
% ADD
% 
% par.pixel_scaling = 0.9992;
% write.subfolder_reco = sprintf('pixscal%.4f',par.pixel_scaling);
% ADD
% 
% par.pixel_scaling = 0.9994;
% write.subfolder_reco = sprintf('pixscal%.4f',par.pixel_scaling);
% ADD
% 
% par.pixel_scaling = 0.9996;
% write.subfolder_reco = sprintf('pixscal%.4f',par.pixel_scaling);
% ADD
% 
% par.pixel_scaling = 0.9998;
% write.subfolder_reco = sprintf('pixscal%.4f',par.pixel_scaling);
% ADD
% 
% par.pixel_scaling = 1.0000;
% write.subfolder_reco = sprintf('pixscal%.4f',par.pixel_scaling);
% ADD
% 
% write.subfolder_reco = '';
% par.pixel_scaling = 0.9985;
% par.ref_range = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p05_reco_loop( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
