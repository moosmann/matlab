function p05_reco_loop_synchroload2017oct_11003440_001( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
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
% Created on 19-Mar-2018 by moosmanj

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

% PREPROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_roi = []; % [y0 y1] vertical roi.  skips first raw_roi(1)-1 lines, reads until raw_roi(2). Not working for *.raw data where images are flipped.
raw_bin = 2; % projection binning factor: 1, 2, or 4
bin_before_filtering = 0; % Binning before pixel filtering is applied, much faster but less effective filtering
excentric_rot_axis = 0; % off-centered rotation axis increasing FOV. -1: left, 0: centeerd, 1: right. influences rot_corr_area1
crop_at_rot_axis = 0; % for recos of scans with excentric rotation axis but WITHOUT projection stitching
stitch_projections = 0; % for 2 pi scans: stitch projection at rotation axis position. Recommended with phase retrieval to reduce artefacts. Standard absorption contrast data should work well without stitching. Subpixel stitching not supported (non-integer rotation axis position is rounded, less/no binning before reconstruction can be used to improve precision).
stitch_method = 'sine'; 'step';'linear'; %  ! CHECK correlation area !
% 'step' : no interpolation, use step function
% 'linear' : linear interpolation of overlap region
% 'sine' : sinusoidal interpolation of overlap region
proj_range = 1; % range of projections to be used (from all found, if empty or 1: all, if scalar: stride, if range: start:incr:end
ref_range = 1; % range of flat fields to be used (from all found), if empty or 1: all. if scalar: stride, if range: start:incr:end
energy = []; % in eV! if empty: read from log file
sample_detector_distance = []; % in m. if empty: read from log file
eff_pixel_size = []; % in m. if empty: read from log file. effective pixel size =  detector pixel size / magnification
dark_FiltPixThresh = [0.01 0.005]; % Dark fields: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
ref_FiltPixThresh = [0.01 0.005]; % Flat fields: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
proj_FiltPixThresh = [0.01 0.005]; % Raw projection: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
correlation_method = 'entropy';'ssim-ml';'diff';'shift';'ssim';'std';'cov';'corr';'cross-entropy12';'cross-entropy21';'cross-entropyx';'none';
% 'ssim-ml' : Matlab's structural similarity index (SSIM), includes Gaussian smoothing
% 'ssim' : own implementation of SSIM, smoothing not yet implemented
% 'entropy' : entropy measure of proj over flat
% 'cov' : cross covariance
% 'corr' : cross correlation = normalized cross covariance
% 'std' : standard deviation of proj over flat
% 'diff': difference of proj and flat
% 'shift': computes relative shift from peak of cross-correlation map
% 'none' : no correlation, use median flat
corr_shift_max_pixelshift = 0.25; % maximum pixelshift allowed for 'shift'-correlation method: if 0 use the best match (i.e. the one with the least shift), if > 0 uses all flats with shifts smaller than corr_shift_max_pixelshift
corr_num_flats = 1; % number of flat fields used for average/median of flats. for 'shift'-correlation its the maximum number
ring_current_normalization = 1; % normalize flat fields and projections by ring current
flat_corr_area1 = [1 floor(100/raw_bin)];%[0.98 1];% correlation area: index vector or relative/absolute position of [first pix, last pix]
flat_corr_area2 = [0.2 0.8]; % correlation area: index vector or relative/absolute position of [first pix, last pix]
decimal_round_precision = 2; % precision when rounding pixel shifts
ring_filter.apply = 0; % ring artifact filter (only for scans without wiggle di wiggle)
ring_filter.apply_before_stitching = 0; % ! Consider when phase retrieval is applied !
ring_filter.method = 'wavelet-fft';'jm';
ring_filter.waveletfft.dec_levels = 2:6; % decomposition levels for 'wavelet-fft'
ring_filter.waveletfft.wname = 'db25';'db30'; % wavelet type for 'wavelet-fft'
ring_filter.waveletfft.sigma = 2.4; %  suppression factor for 'wavelet-fft'
ring_filter.jm.median_width = 11; % [3 11 21 31 39];
% PHASE RETRIEVAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phase_retrieval.apply = 0; % See 'PhaseFilter' for detailed description of parameters !
phase_retrieval.apply_before = 0; % before stitching, interactive mode, etc. For phase-contrast data with an excentric rotation axis phase retrieval should be done afterwards. To find the rotataion axis position use this option in a first run, and then turn it of afterwards.
phase_bin = 1; % Binning factor after phase retrieval, but before tomographic reconstruction
phase_retrieval.method = 'tie';'qpcut';  'tie'; %'qp' 'ctf' 'tie' 'qp2' 'qpcut'
phase_retrieval.reg_par = 2.5; % regularization parameter. larger values tend to blurrier images. smaller values tend to original data.
phase_retrieval.bin_filt = 0.15; % threshold for quasiparticle retrieval 'qp', 'qp2'
phase_retrieval.cutoff_frequ = 1 * pi; % in radian. frequency cutoff in Fourier space for 'qpcut' phase retrieval
phase_retrieval.padding = 1; % padding of intensities before phase retrieval, 0: no padding
% TOMOGRAPHY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
do_tomo = 1; % run tomographic reconstruction
vol_shape = [];%[1.2 1.2 1]; % shape (voxels) of reconstruction volume. in absolute number of voxels or in relative number w.r.t. the default volume which is given by the detector width and height
vol_size = [];%[-1.2 1.2 -1.2 1.2 -1 1]; % 6-component vector [xmin xmax ymin ymax zmin zmax]. if empty, volume is centerd within vol_shape. unit voxel size is assumed. if smaller than 10 values are interpreted as relative size w.r.t. the detector size. Take care bout minus signs!
rot_angle_full = []; % in radians: empty ([]), full angle of rotation, or array of angles. if empty full rotation angles is determined automatically to pi or 2 pi
rot_angle_offset = pi; % global rotation of reconstructed volume
rot_axis_offset = [];% if empty use automatic computation
rot_axis_pos = []; % if empty use automatic computation. either offset or pos has to be empty. can't use both
rot_corr_area1 = []; % ROI to correlate projections at angles 0 & pi. Use [0.75 1] or so for scans with an excentric rotation axis
rot_corr_area2 = []; % ROI to correlate projections at angles 0 & pi
rot_corr_gradient = 0; % use gradient of intensity maps if signal variations are too weak to correlate projections
rot_axis_tilt = 0; % in rad. camera tilt w.r.t rotation axis. if empty calculate from registration of projections at 0 and pi
fbp_filter_type = 'Ram-Lak';'linear'; % see iradonDesignFilter for more options. Ram-Lak according to Kak/Slaney
fpb_filter_freq_cutoff = 1; % Cut-off frequency in Fourier space of the above FBP filter
fbp_filter_padding = 1; % symmetric padding for consistent boundary conditions, 0: no padding
fbp_filter_padding_method = 'symmetric';
butterworth_filter = 0; % use butterworth filter in addition to FBP filter
butterworth_order = 1;
butterworth_cutoff_frequ = 0.8;
astra_pixel_size = 1; % size of a detector pixel: if different from one 'vol_size' needs to be ajusted
take_neg_log = []; % take negative logarithm. if empty, use 1 for attenuation contrast, 0 for phase contrast
% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out_path = '';% absolute path were output data will be stored. !!overwrites the write.to_scratch flag. if empty uses the beamtime directory and either 'processed' or 'scratch_cc'
write.to_scratch = 0; % write to 'scratch_cc' instead of 'processed'
write.flatcor = 0; % save preprocessed flat corrected projections
write.phase_map = 0; % save phase maps (if phase retrieval is not 0)
write.sino = 0; % save sinograms (after preprocessing & before FBP filtering and phase retrieval)
write.phase_sino = 0; % save sinograms of phase maps
write.reco = 1; % save reconstructed slices (if do_tomo=1)
write.float = 1; % single precision (32-bit float) tiff
write.uint16 = 0;
write.uint8 = 0;
reco_bin = 2; % binning factor of reconstructed volume if binned volumes are saved
write.float_binned = 0; % binned single precision (32-bit float) tiff
write.uint16_binned = 0;
write.uint8_binned = 0;
write.uint8_segmented = 0; % experimental: threshold segmentation for histograms with 2 distinct peaks: __/\_/\__
compression_method = 'histo';'full'; 'std'; 'threshold'; % method to compression dynamic range into [0, 1]
compression_parameter = [0.20 0.15]; % compression-method specific parameter
% dynamic range is compressed s.t. new dynamic range assumes
% 'full' : full dynamic range is used
% 'threshold' : [LOW HIGH] = compression_parameter, eg. [-0.01 1]
% 'std' : NUM = compression_parameter, mean +/- NUM*std, dynamic range is rescaled to within -/+ NUM standard deviations around the mean value
% 'histo' : [LOW HIGH] = compression_parameter (100*LOW)% and (100*HIGH)% of the original histogram, e.g. [0.02 0.02]
parfolder = '';% parent folder for 'reco', 'sino', 'phase', and 'flat_corrected'
subfolder_flatcor = ''; % subfolder in 'flat_corrected'
subfolder_phase_map = ''; % subfolder in 'phase_map'
subfolder_sino = ''; % subfolder in 'sino'
subfolder_reco = ''; % subfolder in 'reco'
% INTERACTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
verbose = 1; % print information to standard output
visual_output = 1; % show images and plots during reconstruction
interactive_determination_of_rot_axis = 1; % reconstruct slices with different rotation axis offsets
interactive_determination_of_rot_axis_tilt = 0; % reconstruct slices with different offset AND tilts of the rotation axis
lamino = 0; % find laminography tilt instead camera rotation
fixed_tilt = 0; % fixed other tilt
slice_number = 0.5; % slice number, default: 0.5. if in [0,1): relative, if in (1, N]: absolute
% HARDWARE / SOFTWARE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
use_cluster = 0; % if available: on MAXWELL nodes disp/nova/wga/wgs cluster computation can be used. Recommended only for large data sets since parpool creation and data transfer implies a lot of overhead.
poolsize = 0.60; % number of workers used in a local parallel pool. if 0: use current config. if >= 1: absolute number. if 0 < poolsize < 1: relative amount of all cores to be used. if SLURM scheduling is available, a default number of workers is used.
link_data = 1; % ASTRA data objects become references to Matlab arrays. Reduces memory issues.
gpu_index = []; % GPU Device index to use, Matlab notation: index starts from 1. default: [], uses all
automatic_mode = 0; % Find rotation axis position automatically. NOT IMPLEMENTED!

% Set default. Defines default parameter set to be used.
SET_DEFAULT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETER / DATA SETS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_path = '/asap3/petra3/gpfs/p05/2017/data/11003440/raw/';

scan_path = [raw_path 'syn01_52L_PEEK_12w']; ADD

scan_path = [raw_path 'syn02_52L_PEEK_12w_ccd']; ADD

scan_path = [raw_path 'syn03_52L_PEEK_12w_dcm_cmos_recotest']; ADD

scan_path = [raw_path 'syn04_52L_PEEK_12w_dcm_cmos_recotest']; ADD

scan_path = [raw_path 'syn05_52L_PEEK_12w_dmm_cmos_recotest']; ADD

scan_path = [raw_path 'syn06_52L_PEEK_12w_dmm_ccd_recotest']; ADD

scan_path = [raw_path 'syn07_cor']; ADD

scan_path = [raw_path 'syn08_cor']; ADD

scan_path = [raw_path 'syn09_cor']; ADD

scan_path = [raw_path 'syn10_cor']; ADD

scan_path = [raw_path 'syn11_cor']; ADD

scan_path = [raw_path 'syn12_cor']; ADD

scan_path = [raw_path 'syn13_cor']; ADD

scan_path = [raw_path 'syn15_cor']; ADD

scan_path = [raw_path 'syn16_39L_PEEK_12w']; ADD

scan_path = [raw_path 'syn17_21R_PEEK_8w']; ADD

scan_path = [raw_path 'syn18_88L_Mg10Gd_4w']; ADD

scan_path = [raw_path 'syn19_74R_Mg10Gd_8w']; ADD

scan_path = [raw_path 'syn20_72R_Mg10Gd_12w']; ADD

scan_path = [raw_path 'syn21_82R_Mg10Gd_8w']; ADD

scan_path = [raw_path 'syn22_62L_Mg5Gd_12w']; ADD

scan_path = [raw_path 'syn23_63L_Mg5Gd_12w']; ADD

scan_path = [raw_path 'syn24_78L_Mg5Gd_8w']; ADD

scan_path = [raw_path 'syn25_69R_Mg10Gd_12w']; ADD

scan_path = [raw_path 'syn26_69R_Mg10Gd_12w']; ADD

scan_path = [raw_path 'syn28_69R_Mg10Gd_12w']; ADD

scan_path = [raw_path 'syn29_78L_Mg5Gd_8w']; ADD

scan_path = [raw_path 'syn30_87L_Mg10Gd_4w']; ADD

scan_path = [raw_path 'syn31_97R_Mg10Gd_4w']; ADD

scan_path = [raw_path 'syn32_99R_Mg10Gd_4w']; ADD

scan_path = [raw_path 'syn33_80R_Mg10Gd_8w']; ADD

scan_path = [raw_path 'syn34_79R_Mg10Gd_8w']; ADD

scan_path = [raw_path 'syn35_77R_Mg10Gd_8w']; ADD

scan_path = [raw_path 'syn36_63L_Mg5Gd_12w']; ADD

scan_path = [raw_path 'syn37_69L_Mg10Gd_12w']; ADD

scan_path = [raw_path 'syn38_73R_Mg10Gd_8w']; ADD

scan_path = [raw_path 'syn39_75L_Mg5Gd_8w']; ADD

scan_path = [raw_path 'syn40_69L_Mg10Gd_12w']; ADD

scan_path = [raw_path 'syn41_63L_Mg5Gd_12w']; ADD

scan_path = [raw_path 'syn42_38L_PEEK_8w']; ADD

scan_path = [raw_path 'syn43_38L_PEEK_8w']; ADD

scan_path = [raw_path 'syn44_66L_Mg5Gd_12w']; ADD

scan_path = [raw_path 'syn45_101BL_Mg5Gd_4w']; ADD

scan_path = [raw_path 'syn46_88R_Mg5Gd_4w']; ADD

scan_path = [raw_path 'syn47_100AL_Mg5Gd_4w']; ADD

scan_path = [raw_path 'syn48_89L_Mg5Gd_4w']; ADD

scan_path = [raw_path 'syn49_80L_Mg5Gd_8w']; ADD

scan_path = [raw_path 'syn50_99L_Mg5Gd_4w']; ADD

scan_path = [raw_path 'syn51_87R_Mg5Gd_4w']; ADD

scan_path = [raw_path 'syn52_95R_Mg5Gd_8w']; ADD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3-in-1 sample corrosion, CCD
scan_path = [raw_path 'syn60_cor']; raw_bin = 2; rot_axis_offset = -1.5; ADD

% 3-in-1 sample corrosion, CCD
scan_path = [raw_path 'syn61_cor']; raw_bin = 2; rot_axis_offset = -1.5; ADD

% KIT camera
bin_before_filtering = 1;
SET_DEFAULT;

% Explant, KIT
scan_path = [raw_path 'syn64_91R_Mg5Gd_4w']; raw_bin = 3; rot_axis_offset = 1.75; ADD

scan_path = [raw_path 'syn65_94R_Mg5Gd_8w'];raw_bin = 3; rot_axis_offset = 1.75; ADD

scan_path = [raw_path 'syn66_cor_Punkt2']; raw_bin = 3; rot_axis_offset = 1.75; ADD

scan_path = [raw_path 'syn67_cor_Punkt1_previous_is_Punkt2']; raw_bin = 3; rot_axis_offset = 2.0; ADD

scan_path = [raw_path 'syn68_cor_Punkt3']; raw_bin = 3; rot_axis_offset = 1.83; ADD

scan_path = [raw_path 'syn69_cor_P1_I']; raw_bin = 3; rot_axis_offset = 1.83;  ADD

scan_path = [raw_path 'syn70_cor_P1_2']; raw_bin = 3; rot_axis_offset = 1.83;  ADD

scan_path = [raw_path 'syn71_cor_P1_3']; raw_bin = 3; rot_axis_offset = 2.0; ADD

scan_path = [raw_path 'syn72_cor_P1_4']; raw_bin = 3; rot_axis_offset = 1.83; ADD

scan_path = [raw_path 'syn73_cor_P1_5']; rot_axis_offset = 1.83; ADD

scan_path = [raw_path 'syn74_cor_P1_6']; ADD

scan_path = [raw_path 'syn75_cor_P1_7']; ADD

scan_path = [raw_path 'syn76_cor_P2_1']; ADD

scan_path = [raw_path 'syn77_cor_P2_2']; ADD

scan_path = [raw_path 'syn78_cor_P2_3']; ADD

scan_path = [raw_path 'syn79_cor_P2_4']; ADD

scan_path = [raw_path 'syn80_cor_P2_5']; ADD

scan_path = [raw_path 'syn81_cor_P2_6']; ADD

scan_path = [raw_path 'syn82_cor_P2_7']; ADD

scan_path = [raw_path 'syn83_cor_P3_1']; ADD

scan_path = [raw_path 'syn84_cor_P3_2']; ADD

scan_path = [raw_path 'syn85_cor_P3_3']; ADD

scan_path = [raw_path 'syn86_cor_P3_4']; ADD

scan_path = [raw_path 'syn87_cor_P3_5']; ADD

scan_path = [raw_path 'syn88_cor_P3_6']; ADD

scan_path = [raw_path 'syn89_cor_P3_7']; ADD

scan_path = [raw_path 'syn90_70L']; ADD

scan_path = [raw_path 'syn91_72L']; ADD

scan_path = [raw_path 'syn92_93R']; ADD

scan_path = [raw_path 'syn93_76L']; ADD

scan_path = [raw_path 'syn94_96R']; ADD

scan_path = [raw_path 'syn95_69L']; ADD

scan_path = [raw_path 'syn96_82L']; ADD

scan_path = [raw_path 'syn97_74L']; ADD

scan_path = [raw_path 'syn98_77L']; ADD

scan_path = [raw_path 'syn99_43R']; ADD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p05_reco_loop( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
