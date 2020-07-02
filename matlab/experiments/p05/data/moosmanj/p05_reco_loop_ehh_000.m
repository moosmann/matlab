function p05_reco_loop_script_ehh_000( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
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
% Created on 11-Jan-2018 by moosmanj

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

scan_path = ...
'/asap3/petra3/gpfs/p05/2017/data/11002839/raw/ehh_2017_013_f';
'/asap3/petra3/gpfs/p05/2017/commissioning/c20171218_000_synchro/raw/syn175_cor_Mg5Gd428s_Mg10Gd408s_Mg5Gd8p_30';
'/asap3/petra3/gpfs/p05/2017/commissioning/c20171218_000_synchro/raw/syn174_cor_Mg5Gd416s_Mg10Gd410s_Mg5Gd7p_30';
'/asap3/petra3/gpfs/p05/2017/commissioning/c20171218_000_synchro/raw/syn173_cor_Mg5Gd413s_Mg10Gd409s_Mg5Gd3p_30';
'/asap3/petra3/gpfs/p05/2017/commissioning/c20171218_000_synchro/raw/syn171_cor_Mg5Gd430s_Mg10Gd413s_pmg9p_30';    
'/asap3/petra3/gpfs/p05/2017/commissioning/c20171218_000_synchro/raw/syn170_cor_Mg5Gd414s_Mg10Gd405s_pmg2p_30';
'/asap3/petra3/gpfs/p05/2017/commissioning/c20171218_000_synchro/raw/syn168_cor_Mg5Gd410s_Mg10Gd401s_pmg1p_30';
'/asap3/petra3/gpfs/p05/2017/data/11003527/raw/szeb_86'; % 2 x stitching
'/asap3/petra3/gpfs/p05/2017/data/11003288/raw/syn151_58L_Mg_12_000';
'/asap3/petra3/gpfs/p05/2017/data/11003288/raw/syn156_berit'; % 2 x stitching
'/asap3/petra3/gpfs/p05/2017/data/11003527/raw/szeb_78';
'/asap3/petra3/gpfs/p05/2017/data/11003288/raw/syn166_104R_Mg10Gd_4w_001';
'/asap3/petra3/gpfs/p05/2017/commissioning/c20171115_000_kit_test/raw/fast_test17';
'/asap3/petra3/gpfs/p05/2017/data/11003773/raw/syn132_cor_mg5gd434s_mg10gd408s_12';
'/asap3/petra3/gpfs/p05/2017/data/11003288/raw/syn165_58L_Mg10Gd_12w_000'; %reconstruction is crap
'/asap3/petra3/gpfs/p05/2017/data/11003288/raw/syn166_104R_Mg10Gd_4w_000'; %reconstruction is crap
read_flatcor = 0; % read flatfield-corrected images from disc, skips preprocessing
read_flatcor_path = ''; % subfolder of 'flat_corrected' containing projections
% PREPROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_roi = [800 2200]; % [y0 y1] vertical roi.  skips first raw_roi(1)-1 lines, reads until raw_roi(2). Not working for *.raw data where images are flipped.
raw_bin = 2; % projection binning factor: 1, 2, or 4
bin_before_filtering = 0; % Binning before pixel filtering is applied, much faster but less effective filtering
excentric_rot_axis = 1; % off-centered rotation axis increasing FOV. -1: left, 0: centeerd, 1: right. influences rot_corr_area1
crop_at_rot_axis = 0; % for recos of scans with excentric rotation axis but WITHOUT projection stitching
stitch_projections = 1; % for 2 pi scans: stitch projection at rotation axis position. Recommended with phase retrieval to reduce artefacts. Standard absorption contrast data should work well without stitching. Subpixel stitching not supported (non-integer rotation axis position is rounded, less/no binning before reconstruction can be used to improve precision).
stitch_method = 'linear';'sine'; 'step'; %  ! adjust correlation area if necessary !
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
correlation_method = 'ssim-ml';'entropy';'diff';'shift';'ssim';'std';'cov';'corr';'cross-entropy12';'cross-entropy21';'cross-entropyx';'none';
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
ring_filter = 1; % ring artifact filter (only for scans without wiggle di wiggle)
ring_filter_method = 'jm';'wavelet-fft';
ring_filter_waveletfft_dec_levels = 2:6; % decomposition levels for 'wavelet-fft'
ring_filter_waveletfft_wname = 'db30';'db25'; % wavelet type for 'wavelet-fft'
ring_filter_waveletfft_sigma = 2.4; %  suppression factor for 'wavelet-fft'
ring_filter_jm_median_width = 11; % [3 11 21 31 39];
ringfilt.apply = 1;
ringfilt.method = 'jm';'wavelet-fft';
ringfilt.waveletfft.dec_levels = 2:6; % decomposition levels for 'wavelet-fft'
ringfilter.waveletfft.wname = 'db30';'db25'; % wavelet type for 'wavelet-fft'
ringfilter.waveletfft.sigma = 2.4; %  suppression factor for 'wavelet-fft'
ringfilter.jm.median_width = 11; % [3 11 21 31 39];
% PHASE RETRIEVAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
do_phase_retrieval = 1; % See 'PhaseFilter' for detailed description of parameters !
phase_retrieval_before = 0; % before stitching, interactive mode, etc. For phase-contrast data with an excentric rotation axis phase retrieval should be done afterwards. To find the rotataion axis position use this option in a first run, and then turn it of afterwards.
phase_bin = 2; % Binning factor after phase retrieval, but before tomographic reconstruction
phase_retrieval_method = 'tie';'qpcut';  'tie'; %'qp' 'ctf' 'tie' 'qp2' 'qpcut'
phase_retrieval_reg_par = 2.5; % regularization parameter. larger values tend to blurrier images. smaller values tend to original data.
phase_retrieval_bin_filt = 0.15; % threshold for quasiparticle retrieval 'qp', 'qp2'
phase_retrieval_cutoff_frequ = 1 * pi; % in radian. frequency cutoff in Fourier space for 'qpcut' phase retrieval
phase_padding = 1; % padding of intensities before phase retrieval, 0: no padding
% TOMOGRAPHY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
do_tomo = 1; % run tomographic reconstruction
vol_shape = [];%[1.2 1.2 1]; % shape (voxels) of reconstruction volume. in absolute number of voxels or in relative number w.r.t. the default volume which is given by the detector width and height
vol_size = [];%[-1.2 1.2 -1.2 1.2 -1 1]; % 6-component vector [xmin xmax ymin ymax zmin zmax]. if empty, volume is centerd within vol_shape. unit voxel size is assumed. if smaller than 10 values are interpreted as relative size w.r.t. the detector size. Take care bout minus signs!
rot_angle_full = []; % in radians: empty ([]), full angle of rotation, or array of angles. if empty full rotation angles is determined automatically to pi or 2 pi
rot_angle_offset = pi; % global rotation of reconstructed volume
rot_axis_offset = []; % if empty use automatic computation
rot_axis_pos = []; % if empty use automatic computation. either offset or pos has to be empty. can't use both
rot_corr_area1 = []; % ROI to correlate projections at angles 0 & pi. Use [0.75 1] or so for scans with an excentric rotation axis
rot_corr_area2 = []; % ROI to correlate projections at angles 0 & pi
rot_corr_gradient = 0; % use gradient of intensity maps if signal variations are too weak to correlate projections
rot_axis_tilt = 0; % in rad. camera tilt w.r.t rotation axis. if empty calculate from registration of projections at 0 and pi
fbp_filter_type = 'linear';'Ram-Lak'; % see iradonDesignFilter for more options. Ram-Lak according to Kak/Slaney
fpb_filter_freq_cutoff = 1; % Cut-off frequency in Fourier space of the above FBP filter
fbp_filter_padding = 1; % symmetric padding for consistent boundary conditions, 0: no padding
fbp_filter_padding_method = 'symmetric';
butterworth_filter = 0; % use butterworth filter in addition to FBP filter
butterworth_order = 1;
butterworth_cutoff_frequ = 0.9;
astra_pixel_size = 1; % size of a detector pixel: if different from one 'vol_size' needs to be ajusted
take_neg_log = []; % take negative logarithm. if empty, use 1 for attenuation contrast, 0 for phase contrast
% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out_path = '';% absolute path were output data will be stored. !!overwrites the write_to_scratch flag. if empty uses the beamtime directory and either 'processed' or 'scratch_cc'
write_to_scratch = 1; % write to 'scratch_cc' instead of 'processed'
write_flatcor = 0; % save preprocessed flat corrected projections
write_phase_map = 0; % save phase maps (if phase retrieval is not 0)
write_sino = 0; % save sinograms (after preprocessing & before FBP filtering and phase retrieval)
write_sino_phase = 0; % save sinograms of phase maps
write_reco = 1; % save reconstructed slices (if do_tomo=1)
write_float = 1; % single precision (32-bit float) tiff
write_16bit = 0;
write_8bit = 0;
reco_bin = 2; % binning factor of reconstructed volume if binned volumes are saved
write_float_binned = 0; % binned single precision (32-bit float) tiff
write_16bit_binned = 0;
write_8bit_binned = 0;
write_8bit_segmented = 0; % experimental: threshold segmentation for histograms with 2 distinct peaks: __/\_/\__
compression_method = 'histo';'full'; 'std'; 'threshold'; % method to compression dynamic range into [0, 1]
compression_parameter = [0.20 0.15]; % compression-method specific parameter
% dynamic range is compressed s.t. new dynamic range assumes
% 'full' : full dynamic range is used
% 'threshold' : [LOW HIGH] = compression_parameter, eg. [-0.01 1]
% 'std' : NUM = compression_parameter, mean +/- NUM*std, dynamic range is rescaled to within -/+ NUM standard deviations around the mean value
% 'histo' : [LOW HIGH] = compression_parameter (100*LOW)% and (100*HIGH)% of the original histogram, e.g. [0.02 0.02]
parfolder = '';%sprintf( 'cor_%s', correlation_method);''; % parent folder for 'reco', 'sino', 'phase', and 'flat_corrected'
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
poolsize = 0.60; % number of workers used in a parallel pool. if > 1: absolute number. if 0 < poolsize < 1: relative amount of all cores to be used
link_data = 1; % ASTRA data objects become references to Matlab arrays. Reduces memory issues.
gpu_index = []; % GPU Device index to use, Matlab notation: index starts from 1. default: [], uses all
automatic_mode = 0; % Find rotation axis position automatically. NOT IMPLEMENTED!
automatic_mode_coarse = 'entropy'; %
automatic_mode_fine = 'iso-grad';

% Set default. Defines default parameter set to be used.
SET_DEFAULT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETER / DATA SETS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_path = '/asap3/petra3/gpfs/p05/2017/data/11002839/raw/';

scan_path = [raw_path 'ehh_020_a']; ADD

scan_path = [raw_path 'ehh_020_b']; ADD

scan_path = [raw_path 'ehh_020_c']; ADD

scan_path = [raw_path 'ehh_020_d']; ADD

scan_path = [raw_path 'ehh_020_e']; ADD

scan_path = [raw_path 'ehh_020_f']; ADD

scan_path = [raw_path 'ehh_020_g']; ADD

scan_path = [raw_path 'ehh_020_h']; ADD

scan_path = [raw_path 'ehh_020_i']; ADD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p05_reco_loop( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
