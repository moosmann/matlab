function synchroload2018nov_11005553( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
% Template function to loop over data sets given in the 'PARAMETER / DATA
% SETS' section below. The 'DEFAULT PARAMETERS' section defines the default
% paramters. To add a data / parameter to the loop, define your
% reconstruction parameters followed by 'ADD'. Using 'ADD('r') restores the
% default parameters after a datat set is added. Default parameters are
% defined with 'SET_DEFAULT' (or 'ADD_DEFAULT' or 'ADD('d')' or
% 'ADD('default')'.
%
% Caution: if parameters are changed, they remain changed after a data set
% was added unless 'ADD('r')' ('ADD('restore')) is used which restores 
% all parameters as defined in the 'DEFAULT PARAMETER' section by
% 'SET_DEFAULT'.
%
% Caution: all parameters not defined within the file are taken from main
% reconstruction script 'p05_reco'.
%
% Caution: if default parameters are not defined by closing the DEFAULT
% PARAMETER section with 'SET_DEFAULT', then the first time 'ADD' is called
% defines default parameters.
% 
% Workflow:
% - Copy this file.
% - Check, modify, or copy parameter from 'p05_reco' in  'DEFAULT
%   PARAMETERS' section.
% - Add parameter / data sets to loop over in the ''PARAMETER / DATA SETS'
%   section. 
% - Type 'F5' or call without arguments to list all added data sets.
% - Choose subset of added data sets to loop over by indices (see SUBSETS
%   arguments below).
% - Use RUN_RECO equals 0 with PRINT_PARAMETERS (see below) to check
%   parameters setting.
% - Set RUN_RECO equals 1 to start the loop.
% 
% ARGUMENTS
% SUBSETS : 1D array of integers. subset of added data sets to be looped over
% RUN_RECO : bool. default: 0. 0: loops over the subsets but does not start
% reconstructions, 1: start the reconstruction loop.
% PRINT_PARAMETERS : string or cell of strings, eg {'raw_roi',
% 'parfolder'}. Parameter setting to be printed at each loop iteration.
% Useful in combination with RUN_RECO = 0 to check parameter setting for
% the subset to loop over
%
% Created on 29-Nov-2018 by moosmanj

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

% Define parameter here. Otherwise parameters are taken from the current
% version of 'p05_reco' if not set in the section below.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fast_reco = 0; % !!! OVERWRITES SOME PARAMETERS SET BELOW !!!
stop_after_data_reading(1) = 0; % for data analysis, before flat field correlation 
stop_after_proj_flat_correlation(1) = 0; % for data analysis, after flat field correlation

%%% SCAN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
read_flatcor = 0; % read flatfield-corrected images from disc, skips preprocessing
read_flatcor_path = ''; % subfolder of 'flat_corrected' containing projections
energy = []; % in eV! if empty: read from log file
sample_detector_distance = []; % in m. if empty: read from log file
eff_pixel_size = []; % in m. if empty: read from log file. effective pixel size =  detector pixel size / magnification
pix_scaling = 1;
%%% PREPROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_roi = []; % if []: use full image; if [y0 y1]: vertical ROI, skips first raw_roi(1)-1 lines, reads until raw_roi(2). When raw_roi(2) < 0 reads until end - |raw_roi(2)|; if negative scalar: auto roi, selects ROI automatically.Not working for *.raw data where images are flipped.
raw_bin = 2; % projection binning factor: integer
im_trafo = ''; % string to be evaluated after reading data in the case the image is flipped/rotated/etc due to changes at the beamline, e.g. 'rot90(im)'
bin_before_filtering(1) = 0; % Apply binning before filtering pixel. less effective, but much faster especially for KIT camera.
excentric_rot_axis = 0; % off-centered rotation axis increasing FOV. -1: left, 0: centeerd, 1: right. influences tomo.rot_axis.corr_area1
crop_at_rot_axis = 0; % for recos of scans with excentric rotation axis but WITHOUT projection stitching
stitch_projections = 0; % for 2 pi scans: stitch projection at rotation axis position. Recommended with phase retrieval to reduce artefacts. Standard absorption contrast data should work well without stitching. Subpixel stitching not supported (non-integer rotation axis position is rounded, less/no binning before reconstruction can be used to improve precision).
stitch_method = 'sine'; 'step';'linear'; %  ! CHECK correlation area !
% 'step' : no interpolation, use step function
% 'linear' : linear interpolation of overlap region
% 'sine' : sinusoidal interpolation of overlap region
proj_range = []; % range of projections to be used (from all found, if empty or 1: all, if scalar: stride, if range: start:incr:end
ref_range = []; % range of flat fields to be used (from all found), if empty or 1: all. if scalar: stride, if range: start:incr:end
pixel_filter_threshold_dark = [0.01 0.005]; % Dark fields: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_threshold_flat = [0.01 0.005]; % Flat fields: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_threshold_proj = [0.01 0.005]; % Raw projection: threshold parameter for hot/dark pixel filter, for details see 'FilterPixe
ring_current_normalization = 1; % normalize flat fields and projections by ring current
image_correlation.method = 'ssim-ml';'entropy';'diff';'shift';'ssim';'std';'cov';'corr';'cross-entropy12';'cross-entropy21';'cross-entropyx';
% Correlation of projections and flat fields. Essential for DCM data. Even
% though less efficient for DMM data, it usually improves reconstruction quality.
% Available methods ('ssim-ml'/'entropy' usually work best):
% 'none' : no correlation, uses median flat, for fast recos
% 'ssim-ml' : Matlab's structural similarity index (SSIM), includes Gaussian smoothing
% 'entropy' : entropy measure of proj over flat
% 'ssim' : own implementation of SSIM, smoothing not yet implemented, usually worse, but sometimes better than 'ssim-ml'
% 'cov' : cross covariance
% 'corr' : cross correlation = normalized cross covariance
% 'std' : standard deviation of proj over flat
% 'diff': difference of proj and flat
% 'shift': computes relative shift from peak of cross-correlation map
% 'cross-entropy*' : variants of (asymmetric) cross entropy
image_correlation.num_flats = 1; % number of flat fields used for average/median of flats. for 'shift'-correlation its the maximum number
image_correlation.area_width = [0 0.02];%[0.98 1];% correlation area: index vector or relative/absolute position of [first pix, last pix]
image_correlation.area_height = [0.2 0.8]; % correlation area: index vector or relative/absolute position of [first pix, last pix]
image_correlation.shift.max_pixelshift = 0.25; % maximum pixelshift allowed for 'shift'-correlation method: if 0 use the best match (i.e. the one with the least shift), if > 0 uses all flats with shifts smaller than image_correlation.shift.max_pixelshift
ring_filter.apply = 0; % ring artifact filter (use only for scans without lateral sample movement)
ring_filter.apply_before_stitching = 0; % ! Consider when phase retrieval is applied !
ring_filter.method = 'jm'; 'wavelet-fft';
ring_filter.waveletfft.dec_levels = 2:5; % decomposition levels for 'wavelet-fft'
ring_filter.waveletfft.wname = 'db25';'db30'; % wavelet type, see 'FilterStripesCombinedWaveletFFT' or 'waveinfo'
ring_filter.waveletfft.sigma = 2.4; %  suppression factor for 'wavelet-fft'
ring_filter.jm.median_width = 11; % multiple widths are applied consecutively, eg [3 11 21 31 39];
strong_abs_thresh = 1; % if 1: does nothing, if < 1: flat-corrected values below threshold are set to one. Try with algebratic reco techniques.
decimal_round_precision = 2; % precision when rounding pixel shifts
%%% PHASE RETRIEVAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phase_retrieval.apply = 0; % See 'PhaseFilter' for detailed description of parameters !
phase_retrieval.apply_before = 0; % before stitching, interactive mode, etc. For phase-contrast data with an excentric rotation axis phase retrieval should be done afterwards. To find the rotataion axis position use this option in a first run, and then turn it of afterwards.
phase_retrieval.post_binning_factor = 1; % Binning factor after phase retrieval, but before tomographic reconstruction
phase_retrieval.method = 'tie';'qpcut'; %'qp' 'ctf' 'tie' 'qp2' 'qpcut'
phase_retrieval.reg_par = 1.5; % regularization parameter. larger values tend to blurrier images. smaller values tend to original data.
phase_retrieval.bin_filt = 0.15; % threshold for quasiparticle retrieval 'qp', 'qp2'
phase_retrieval.cutoff_frequ = 2 * pi; % in radian. frequency cutoff in Fourier space for 'qpcut' phase retrieval
phase_retrieval.padding = 1; % padding of intensities before phase retrieval, 0: no padding
%%% TOMOGRAPHY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tomo.run = 1; % run tomographic reconstruction
tomo.run_interactive_mode = 1; % if tomo.run = 0, intendet to determine rot axis positions without processing the full tomogram;
tomo.reco_mode = '3D'; 'slice'; % slice-wise or full 3D backprojection. 'slice': volume must be centered at origin & no support of rotation axis tilt, reco binning, save compressed
tomo.vol_size = []; %[-0.5 0.5 -0.5 0.5 -0.5 0.5];% 6-component vector [xmin xmax ymin ymax zmin zmax], for excentric rot axis pos / extended FoV;. if empty, volume is centerd within tomo.vol_shape. unit voxel size is assumed. if smaller than 10 values are interpreted as relative size w.r.t. the detector size. Take care bout minus signs! Note that if empty vol_size is dependent on the rotation axis position.
tomo.vol_shape = []; %[1 1 1] shape (# voxels) of reconstruction volume. used for excentric rot axis pos. if empty, inferred from 'tomo.vol_size'. in absolute numbers of voxels or in relative number w.r.t. the default volume which is given by the detector width and height.
tomo.rot_angle.full_range = []; % in radians: empty ([]), full angle of rotation, or array of angles. if empty full rotation angles is determined automatically to pi or 2 pi
tomo.rot_angle.offset = pi; % global rotation of reconstructed volume
tomo.rot_axis.offset = []; % if empty use automatic computation
tomo.rot_axis.offset_shift_range = []; %[]; % absolute lateral movement in pixels during fly-shift-scan, overwrite lateral shift read out from hdf5 log
tomo.rot_axis.position = []; % if empty use automatic computation. EITHER OFFSET OR POSITION MUST BE EMPTY. YOU MUST NOT USE BOTH!
tomo.rot_axis.tilt = 0; % in rad. camera tilt w.r.t rotation axis. if empty calculate from registration of projections at 0 and pi
tomo.rot_axis.corr_area1 = []; % ROI to correlate projections at angles 0 & pi. Use [0.75 1] or so for scans with an excentric rotation axis
tomo.rot_axis.corr_area2 = []; % ROI to correlate projections at angles 0 & pi
tomo.rot_axis.corr_gradient = 0; % use gradient of intensity maps if signal variations are too weak to correlate projections
tomo.fbp_filter.type = 'Ram-Lak';'linear'; % see iradonDesignFilter for more options. Ram-Lak according to Kak/Slaney
tomo.fbp_filter.freq_cutoff = 1; % Cut-off frequency in Fourier space of the above FBP filter
tomo.fbp_filter.padding = 1; % symmetric padding for consistent boundary conditions, 0: no padding
tomo.fbp_filter.padding_method = 'symmetric';
tomo.butterworth_filter.apply = 0; % use butterworth filter in addition to FBP filter
tomo.butterworth_filter.order = 1;
tomo.butterworth_filter.frequ_cutoff = 0.9;
tomo.astra_pixel_size = 1; % detector pixel size for reconstruction: if different from one 'tomo.vol_size' must to be ajusted, too!
tomo.take_neg_log = []; % take negative logarithm. if empty, use 1 for attenuation contrast, 0 for phase contrast
tomo.algorithm = 'fbp';'sirt'; 'cgls';'sart';'em';'fbp-astra'; % SART/EM only work for 3D reco mode
tomo.iterations = 40; % for 'sirt' or 'cgls'.
tomo.sirt.MinConstraint = []; % If specified, all values below MinConstraint will be set to MinConstraint. This can be used to enforce non-negative reconstructions, for example.
tomo.sirt.MaxConstraint = []; % If specified, all values above MaxConstraint will be set to MaxConstraint.
%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write.path = '';'/gpfs/petra3/scratch/moosmanj';% absolute path were output data will be stored. !!overwrites the write.to_scratch flag. if empty uses the beamtime directory and either 'processed' or 'scratch_cc'
write.to_scratch = 0; % write to 'scratch_cc' instead of 'processed'
write.parfolder = '';% parent folder to 'reco', 'sino', 'phase', and 'flat_corrected'
write.subfolder.flatcor = ''; % subfolder in 'flat_corrected'
write.subfolder.phase_map = ''; % subfolder in 'phase_map'
write.subfolder.sino = ''; % subfolder in 'sino'
write.subfolder.reco = ''; % subfolder in 'reco'
write.flatcor = 1; % save preprocessed flat corrected projections
write.flatcor_shift_cropped = 1; % save lateral shift corrected projections, projections are not interpolated, but cropped to nearest integer pixel
write.phase_map = 0; % save phase maps (if phase retrieval is not 0)
write.sino = 1; % save sinograms (after preprocessing & before FBP filtering and phase retrieval)
write.sino_shift_cropped = 1; % save cropped sinos without lateral shift
write.phase_sino = 0; % save sinograms of phase maps
write.reco = 1; % save reconstructed slices (if tomo.run=1)
write.float = 1; % single precision (32-bit float) tiff
write.uint16 = 0; % save 16bit unsigned integer tiff using 'write.compression.method'
write.uint8 = 0; % save binned 8bit unsigned integer tiff using 'write.compression.method'
% Optionally save binned reconstructions, only works in '3D' reco_mode
write.float_binned = 1; % save binned single precision (32-bit float) tiff
write.uint16_binned = 0; % save binned 16bit unsigned integer tiff using 'write.compression.method'
write.uint8_binned = 0; % save binned 8bit unsigned integer tiff using 'wwrite.compression.method'
write.reco_binning_factor = 2; % IF BINNED VOLUMES ARE SAVED: binning factor of reconstructed volume
write.compression.method =  'outlier';'threshold';'histo';'full'; 'std'; % method to compression dynamic range into [0, 1]
write.compression.parameter = [0.02 0.02]; % compression-method specific parameter
% dynamic range is compressed s.t. new dynamic range assumes
% 'outlier' : [LOW, HIGH] = write.compression.parameter, eg. [0.01 0.03], outlier given in percent, if scalear LOW = HIGH.
% 'full' : full dynamic range is used
% 'threshold' : [LOW HIGH] = write.compression.parameter, eg. [-0.01 1]
% 'std' : NUM = write.compression.parameter, mean +/- NUM*std, dynamic range is rescaled to within -/+ NUM standard deviations around the mean value
% 'histo' : [LOW HIGH] = write.compression.parameter (100*LOW)% and (100*HIGH)% of the original histogram, e.g. [0.02 0.02]
%%% INTERACTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
verbose = 1; % print information to standard output
visual_output = 1; % show images and plots during reconstruction
interactive_mode.rot_axis_pos = 1; % reconstruct slices with dif+ferent rotation axis offsets
interactive_mode.rot_axis_tilt = 0; % reconstruct slices with different offset AND tilts of the rotation axis
interactive_mode.lamino = 0; % find laminography tilt instead camera rotation
interactive_mode.fixed_other_tilt = 0; % fixed other tilt
interactive_mode.slice_number = 0.5; % default slice number. if in [0,1): relative, if in (1, N]: absolute
interactive_mode.phase_retrieval = 0; % Interactive retrieval to determine regularization parameter
%%% HARDWARE / SOFTWARE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
use_cluster = 0; % if available: on MAXWELL nodes disp/nova/wga/wgs cluster computation can be used. Recommended only for large data sets since parpool creation and data transfer implies a lot of overhead.
poolsize = 0.60; % number of workers used in a local parallel pool. if 0: use current config. if >= 1: absolute number. if 0 < poolsize < 1: relative amount of all cores to be used. if SLURM scheduling is available, a default number of workers is used.
tomo.astra_link_data = 1; % ASTRA data objects become references to Matlab arrays. Reduces memory issues.
tomo.astra_gpu_index = []; % GPU Device index to use, Matlab notation: index starts from 1. default: [], uses all
%%% EXPERIMENTAL OR NOT YET IMPLEMENTED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write.uint8_segmented = 0; % experimental: threshold segmentation for histograms with 2 distinct peaks: __/\_/\__

% Set default. Defines default parameter set to be used.
SET_DEFAULT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETER / DATA SETS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_path = '/asap3/petra3/gpfs/p05/2018/data/11005553/raw/';

scan_path = [raw_path 'syn002_64R_Mg10Gd_12w']; ADD
scan_path = [raw_path 'syn003_64R_Mg10Gd_12w_kit']; ADD

%% Load, CCD
raw_roi = 2;
scan_path = [raw_path 'syn004_24R_PEEK_8w_000_radio']; ADD
% strong movement
tomo.rot_axis.offset = -0.25 * 2 / raw_bin;
% beam current drops to 80 mA
tomo.rot_axis.offset = -0.1 * 2 / raw_bin;
scan_path = [raw_path 'syn004_24R_PEEK_8w_000']; ADD
% small local movement artefacts
tomo.rot_axis.offset = 0.1 * 2 / raw_bin;
scan_path = [raw_path 'syn004_24R_PEEK_8w_001']; ADD
tomo.rot_axis.offset = 0.1 * 2 / raw_bin;
scan_path = [raw_path 'syn004_24R_PEEK_8w_002']; ADD
tomo.rot_axis.offset = 0.0 * 2 / raw_bin;
scan_path = [raw_path 'syn004_24R_PEEK_8w_003']; ADD
tomo.rot_axis.offset = 0.1 * 2 / raw_bin;
scan_path = [raw_path 'syn004_24R_PEEK_8w_004']; ADD
tomo.rot_axis.offset = 0.0 * 2 / raw_bin;
scan_path = [raw_path 'syn004_24R_PEEK_8w_005']; ADD


scan_path = [raw_path 'syn004_24R_PEEK_8w_006']; ADD

scan_path = [raw_path 'syn005_94R_Mg5Gd_8w']; ADD

scan_path = [raw_path 'syn007_24R_PEEK_8w']; ADD
scan_path = [raw_path 'syn008_24R_PEEK_8w_47N']; ADD
scan_path = [raw_path 'syn009_24R_PEEK_8w_crashed_relaxed_backto30N']; ADD

scan_path = [raw_path 'syn010_34L_PEEK_8w']; ADD

%% Load, CCD
scan_path = [raw_path 'syn011_105L_Mg5Gd_4w_000_radio']; ADD
scan_path = [raw_path 'syn011_105L_Mg5Gd_4w_000']; ADD
scan_path = [raw_path 'syn012_105L_Mg5Gd_4w_000']; ADD
tomo.rot_axis.offset = -7 * 2 / raw_bin;
scan_path = [raw_path 'syn013_105L_Mg5Gd_4w_000']; ADD
scan_path = [raw_path 'syn013_105L_Mg5Gd_4w_001']; ADD
scan_path = [raw_path 'syn013_105L_Mg5Gd_4w_002']; ADD
scan_path = [raw_path 'syn013_105L_Mg5Gd_4w_003']; ADD
scan_path = [raw_path 'syn013_105L_Mg5Gd_4w_004']; ADD
scan_path = [raw_path 'syn013_105L_Mg5Gd_4w_005']; ADD
scan_path = [raw_path 'syn013_105L_Mg5Gd_4w_006']; ADD
scan_path = [raw_path 'syn013_105L_Mg5Gd_4w_007']; ADD
scan_path = [raw_path 'syn013_105L_Mg5Gd_4w_008']; ADD
scan_path = [raw_path 'syn013_105L_Mg5Gd_4w_009']; ADD

interactive_mode.rot_axis_pos = 1;
tomo.rot_axis.offset = 1 * 2 / raw_bin;
scan_path = [raw_path 'syn014_105R_Mg10Gd_4w_000']; ADD
tomo.rot_axis.offset = -1.1 * 2 / raw_bin;
scan_path = [raw_path 'syn014_105R_Mg10Gd_4w_001']; ADD
tomo.rot_axis.offset = 0.2 * 2 / raw_bin;
scan_path = [raw_path 'syn014_105R_Mg10Gd_4w_002']; ADD
tomo.rot_axis.offset = 0.1 * 2 / raw_bin;
scan_path = [raw_path 'syn014_105R_Mg10Gd_4w_003']; ADD
tomo.rot_axis.offset = 0.1 * 2 / raw_bin;
scan_path = [raw_path 'syn014_105R_Mg10Gd_4w_004']; ADD
tomo.rot_axis.offset = 0.1 * 2 / raw_bin;
scan_path = [raw_path 'syn014_105R_Mg10Gd_4w_005']; ADD
tomo.rot_axis.offset = 0.0 * 2 / raw_bin;
scan_path = [raw_path 'syn014_105R_Mg10Gd_4w_006']; ADD
tomo.rot_axis.offset = 0.2 * 2 / raw_bin;
scan_path = [raw_path 'syn014_105R_Mg10Gd_4w_007']; ADD
tomo.rot_axis.offset = 0.6 * 2 / raw_bin;
scan_path = [raw_path 'syn014_105R_Mg10Gd_4w_008']; ADD
% small locally varying movement artefacts
tomo.rot_axis.offset = 0.6 * 2 / raw_bin;
scan_path = [raw_path 'syn014_105R_Mg10Gd_4w_009']; ADD

% good quality
interactive_mode.rot_axis_pos = 0;
raw_bin = 3;
tomo.rot_axis.offset = 0.05 * 2 / raw_bin;
scan_path = [raw_path 'syn015_86L_Mg10Gd_4w_beforedrilling']; ADD

scan_path = [raw_path 'syn016_67L_Mg10Gd_8w_000']; ADD

% sligtly worse than before drilling
raw_bin = 3;
tomo.rot_axis.offset = 0.05 * 2 / raw_bin;
scan_path = [raw_path 'syn017_86L_Mg10Gd_4w_afterdrilling']; ADD

% strong movement artefacts
tomo.rot_axis.offset = -0.95 * 2 / raw_bin; % metal pin
%tomo.rot_axis.offset = 1.25 * 2 / raw_bin; % sample holder bottom
scan_path = [raw_path 'syn018_86L_Mg10Gd_4w_afterdrilling_push_000']; ADD
% movement artefacts
scan_path = [raw_path 'syn018_86L_Mg10Gd_4w_afterdrilling_push_001']; ADD
tomo.rot_axis.offset = -0.75 * 2 / raw_bin;
scan_path = [raw_path 'syn018_86L_Mg10Gd_4w_afterdrilling_push_002']; ADD
% local movement artefacts
tomo.rot_axis.offset = -0.15 * 2 / raw_bin;
scan_path = [raw_path 'syn018_86L_Mg10Gd_4w_afterdrilling_push_003']; ADD
% local movement artefacts
tomo.rot_axis.offset = -0.025 * 2 / raw_bin;
scan_path = [raw_path 'syn018_86L_Mg10Gd_4w_afterdrilling_push_004']; ADD

% good
tomo.rot_axis.offset = -0.45 * 2 / raw_bin;
scan_path = [raw_path 'syn019_86L_Mg10Gd_4w_afterpushing']; ADD


scan_path = [raw_path 'syn020_LongTerm_18003_bot']; ADD
scan_path = [raw_path 'syn020_LongTerm_18003_mid']; ADD
scan_path = [raw_path 'syn020_LongTerm_18003_top']; ADD

scan_path = [raw_path 'syn021_LongTerm_18005_bot']; ADD
scan_path = [raw_path 'syn021_LongTerm_18005_mid']; ADD
scan_path = [raw_path 'syn021_LongTerm_18005_top']; ADD

scan_path = [raw_path 'syn022_LongTerm_18008_bot']; ADD
scan_path = [raw_path 'syn022_LongTerm_18008_mid']; ADD
scan_path = [raw_path 'syn022_LongTerm_18008_top']; ADD

scan_path = [raw_path 'syn023_13R_PEEK_4w']; ADD

%% To do
interactive_mode.rot_axis_pos = 1;
% 0-180: 0.75
scan_path = [raw_path 'syn026_femur_55L_000']; ADD
scan_path = [raw_path 'syn026_femur_55L_001']; ADD

% good, some small local movement
tomo.rot_axis.offset = 1.2 * 3 / raw_bin;
scan_path = [raw_path 'syn027_20R_PEEK_4w']; ADD

% quite some movment in on half of the sample
tomo.rot_axis.offset = 1.5 * 3 / raw_bin;
scan_path = [raw_path 'syn028_60L_Mg10Gd_12w']; ADD

% Load
scan_path = [raw_path 'syn029_84R_Mg5Gd_4w_000']; ADD
scan_path = [raw_path 'syn029_84R_Mg5Gd_4w_001']; ADD
scan_path = [raw_path 'syn029_84R_Mg5Gd_4w_002']; ADD
scan_path = [raw_path 'syn029_84R_Mg5Gd_4w_003']; ADD
% technical problems, scan continues
scan_path = [raw_path 'syn029_84R_Mg5Gd_4w_restart_004']; ADD
scan_path = [raw_path 'syn029_84R_Mg5Gd_4w_restart_005']; ADD
scan_path = [raw_path 'syn029_84R_Mg5Gd_4w_restart_006']; ADD
scan_path = [raw_path 'syn030_84R_Mg5Gd_4w_restart_pushed']; ADD
scan_path = [raw_path 'syn032_84R_Mg5Gd_4w_restart_pushed']; ADD

% Good quality
raw_bin = 3; % CMOS
tomo.rot_axis.offset = 2.65 * 3 / raw_bin;
scan_path = [raw_path 'syn033_68R_Mg10Gd_12w']; ADD

scan_path = [raw_path 'syn034_59R_Mg5Gd_12w']; ADD

scan_path = [raw_path 'syn035_56L_Mg10Gd_12w']; ADD

scan_path = [raw_path 'syn036_LongTerm_18007_bottom']; ADD
scan_path = [raw_path 'syn036_LongTerm_18007_middle']; ADD
scan_path = [raw_path 'syn036_LongTerm_18007_top']; ADD

scan_path = [raw_path 'syn037_LongTerm_18006_bottom']; ADD
scan_path = [raw_path 'syn037_LongTerm_18006_middle']; ADD
scan_path = [raw_path 'syn037_LongTerm_18006_top']; ADD

scan_path = [raw_path 'syn038_LongTerm_18002_bottom']; ADD
scan_path = [raw_path 'syn038_LongTerm_18002_middle']; ADD
scan_path = [raw_path 'syn038_LongTerm_18002_top']; ADD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p05_reco_loop( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
