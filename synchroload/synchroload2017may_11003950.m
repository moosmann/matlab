function synchroload2017may_11003950( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
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
fast_reco = 0; % !!! OVERWRITES SOME PARAMETERS SET BELOW !!!
stop_after_data_reading(1) = 0; % for data analysis, before flat field correlation
stop_after_proj_flat_correlation(1) = 0; % for data analysis, after flat field correlation

%%% SCAN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
read_flatcor = 0; % read preprocessed flatfield-corrected projections. CHECK if negative log has to be taken!
read_flatcor_path = ''; % subfolder of 'flat_corrected' containing projections
read_sino = 0; % read preprocessed sinograms. CHECK if negative log has to be taken!
read_sino_folder = 'trans04'; % subfolder to scan path
energy = []; % in eV! if empty: read from log file
sample_detector_distance = []; % in m. if empty: read from log file
eff_pixel_size = []; % in m. if empty: read from log file. effective pixel size =  detector pixel size / magnification
pix_scaling = 1; % to account for beam divergence if pixel size was determined (via MTF) at single distance only
%%% PREPROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_roi = -8; % vertical and/or horizontal ROI; (1,1) coordinate = top left pixel; supports absolute, relative, negative, and mixed indexing.
raw_bin = 2; % projection binning factor: integer
im_trafo = ''; % string to be evaluated after reading data in the case the image is flipped/rotated/etc due to changes at the beamline, e.g. 'rot90(im)'
bin_before_filtering(1) = 0; % Apply binning before filtering pixel. less effective, but much faster especially for KIT camera.
excentric_rot_axis = 0; % off-centered rotation axis increasing FOV. -1: left, 0: centeerd, 1: right. influences tomo.rot_axis.corr_area1
crop_at_rot_axis = 0; % for recos of scans with excentric rotation axis but WITHOUT projection stitching
stitch_projections = 0; % for 2 pi scans: stitch projection at rotation axis position. Recommended with phase retrieval to reduce artefacts. Standard absorption contrast data should work well without stitching. Subpixel stitching not supported (non-integer rotation axis position is rounded, less/no binning before reconstruction can be used to improve precision).
stitch_method = 'sine'; 'step';'linear'; %  ! CHECK correlation area !
proj_range = []; % range of projections to be used (from all found, if empty or 1: all, if scalar: stride, if range: start:incr:end
ref_range = []; % range of flat fields to be used (from all found), if empty or 1: all. if scalar: stride, if range: start:incr:end
pixel_filter_threshold_dark = [0.01 0.005]; % Dark fields: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_threshold_flat = [0.01 0.005]; % Flat fields: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_threshold_proj = [0.01 0.005]; % Raw projection: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_radius = []; % Default: [3 3]. Increase only if blobs of zeros or other artefacts are expected. Can increase processing time heavily.
ring_current_normalization = 1; % normalize flat fields and projections by ring current
image_correlation.method = 'ssim-ml';'none';'entropy';'entropy';'none';'diff';'shift';'ssim';'std';'cov';'corr';'cross-entropy12';'cross-entropy21';'cross-entropyx';
image_correlation.num_flats = 1; % number of flat fields used for average/median of flats. for 'shift'-correlation its the maximum number
image_correlation.area_width = [0 0.01];%[0.98 1];% correlation area: index vector or relative/absolute position of [first pix, last pix]
image_correlation.area_height = [0.3 0.7]; % correlation area: index vector or relative/absolute position of [first pix, last pix]
image_correlation.shift.max_pixelshift = 0.25; % maximum pixelshift allowed for 'shift'-correlation method: if 0 use the best match (i.e. the one with the least shift), if > 0 uses all flats with shifts smaller than image_correlation.shift.max_pixelshift
ring_filter.apply = 1; % ring artifact filter (use only for scans without lateral sample movement)
ring_filter.apply_before_stitching = 0; % ! Consider when phase retrieval is applied !
ring_filter.method = 'jm'; 'wavelet-fft';
ring_filter.waveletfft.dec_levels = 7; % decomposition levels for 'wavelet-fft'
ring_filter.waveletfft.wname = 'db25';'db30'; % wavelet type, see 'FilterStripesCombinedWaveletFFT' or 'waveinfo'
ring_filter.waveletfft.sigma = 2.4; %  suppression factor for 'wavelet-fft'
ring_filter.jm.median_width = 11; % multiple widths are applied consecutively, eg [3 11 21 31 39];
strong_abs_thresh = 1; % if 1: does nothing, if < 1: flat-corrected values below threshold are set to one. Try with algebratic reco techniques.
decimal_round_precision = 2; % precision when rounding pixel shifts
%%% PHASE RETRIEVAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phase_retrieval.apply = 0; % See 'PhaseFilter' for detailed description of parameters !
phase_retrieval.apply_before = 1; % before stitching, interactive mode, etc. For phase-contrast data with an excentric rotation axis phase retrieval should be done afterwards. To find the rotataion axis position use this option in a first run, and then turn it of afterwards.
phase_retrieval.post_binning_factor = 1; % Binning factor after phase retrieval, but before tomographic reconstruction
phase_retrieval.method = 'tie';'qp';'qpcut'; %'qp' 'ctf' 'tie' 'qp2' 'qpcut'
phase_retrieval.reg_par = 1.5; % regularization parameter. larger values tend to blurrier images. smaller values tend to original data.
phase_retrieval.bin_filt = 0.1; % threshold for quasiparticle retrieval 'qp', 'qp2'
phase_retrieval.cutoff_frequ = 2 * pi; % in radian. frequency cutoff in Fourier space for 'qpcut' phase retrieval
phase_retrieval.padding = 1; % padding of intensities before phase retrieval, 0: no padding
%%% TOMOGRAPHY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tomo.run = 1; % run tomographic reconstruction
tomo.run_interactive_mode = 1; % if tomo.run = 0, intendet to determine rot axis positions without processing the full tomogram;
tomo.reco_mode = '3D'; 'slice'; % slice-wise or full 3D backprojection. 'slice': volume must be centered at origin & no support of rotation axis tilt, reco binning, save compressed
tomo.vol_size = [];%[-1 1 -1 1 -0.5 0.5];% 6-component vector [xmin xmax ymin ymax zmin zmax], for excentric rot axis pos / extended FoV;. if empty, volume is centerd within tomo.vol_shape. unit voxel size is assumed. if smaller than 10 values are interpreted as relative size w.r.t. the detector size. Take care bout minus signs! Note that if empty vol_size is dependent on the rotation axis position.
tomo.vol_shape = []; %[1 1 1] shape (# voxels) of reconstruction volume. used for excentric rot axis pos. if empty, inferred from 'tomo.vol_size'. in absolute numbers of voxels or in relative number w.r.t. the default volume which is given by the detector width and height.
tomo.rot_angle.full_range = 2*pi; % in radians. if []: full angle of rotation including additional increment, or array of angles. if empty full rotation angles is determined automatically to pi or 2 pi
tomo.rot_angle.offset = pi; % global rotation of reconstructed volume
tomo.rot_axis.offset = 0.7;%[]; % if empty use automatic computation
tomo.rot_axis.offset_shift_range = []; %[]; % absolute lateral movement in pixels during fly-shift-scan, overwrite lateral shift read out from hdf5 log
tomo.rot_axis.position = []; % if empty use automatic computation. EITHER OFFSET OR POSITION MUST BE EMPTY. YOU MUST NOT USE BOTH!
tomo.rot_axis.tilt = 0.001; % in rad. camera tilt w.r.t rotation axis. if empty calculate from registration of projections at 0 and pi
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
write.path = ''; %'/gpfs/petra3/scratch/moosmanj'; % absolute path were output data will be stored. !!overwrites the write.to_scratch flag. if empty uses the beamtime directory and either 'processed' or 'scratch_cc'
write.to_scratch = 0; % write to 'scratch_cc' instead of 'processed'
write.parfolder = '';% parent folder to 'reco', 'sino', 'phase', and 'flat_corrected'
write.subfolder.flatcor = ''; % subfolder in 'flat_corrected'
write.subfolder.phase_map = ''; % subfolder in 'phase_map'
write.subfolder.sino = ''; % subfolder in 'sino'
write.subfolder.reco = ''; % subfolder in 'reco'
write.flatcor = 1; % save preprocessed flat corrected projections
write.flatcor_shift_cropped = 1; % save lateral shift corrected projections, projections are not interpolated, but cropped to nearest integer pixel
write.phase_map = 0; % save phase maps (if phase retrieval is not 0)
write.sino = 0; % save sinograms (after preprocessing & before FBP filtering and phase retrieval)
write.sino_shift_cropped = 0; % save cropped sinos without lateral shift
write.phase_sino = 0; % save sinograms of phase maps
write.reco = 1; % save reconstructed slices (if tomo.run=1)
write.float = 1; % single precision (32-bit float) tiff
write.uint16 = 0; % save 16bit unsigned integer tiff using 'write.compression.method'
write.uint8 = 0; % save binned 8bit unsigned integer tiff using 'write.compression.method'
% Optionally save binned reconstructions, only works in '3D' reco_mode
write.float_binned = 0; % save binned single precision (32-bit float) tiff
write.uint16_binned = 0; % save binned 16bit unsigned integer tiff using 'write.compression.method'
write.uint8_binned = 0; % save binned 8bit unsigned integer tiff using 'wwrite.compression.method'
write.reco_binning_factor = 2; % IF BINNED VOLUMES ARE SAVED: binning factor of reconstructed volume
write.compression.method = 'outlier';'histo';'full'; 'std'; 'threshold'; % method to compression dynamic range into [0, 1]
write.compression.parameter = [0.02 0.02]; % compression-method specific parameter
%%% INTERACTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
verbose = 1; % print information to standard output
visual_output = 1; % show images and plots during reconstruction
interactive_mode.rot_axis_pos = 1; % reconstruct slices with dif+ferent rotation axis offsets
interactive_mode.rot_axis_tilt = 0; % reconstruct slices with different offset AND tilts of the rotation axis
interactive_mode.lamino = 0; % find laminography tilt instead camera rotation
interactive_mode.fixed_other_tilt = 0; % fixed other tilt
interactive_mode.slice_number = 0.5; % default slice number. if in [0,1): relative, if in (1, N]: absolute
interactive_mode.phase_retrieval = 0; % I
use_cluster = 0; % if available: on MAXWELL nodes disp/nova/wga/wgs cluster computation can be used. Recommended only for large data sets since parpool creation and data transfer implies a lot of overhead.
poolsize = 0.7; % number of workers used in a local parallel pool. if 0: use current config. if >= 1: absolute number. if 0 < poolsize < 1: relative amount of all cores to be used. if SLURM scheduling is available, a default number of workers is used.
tomo.astra_link_data = 1; % ASTRA data objects become references to Matlab arrays. Reduces memory issues.
tomo.astra_gpu_index = []; % GPU Device index to use, Matlab notation: index starts from 1. default: [], uses all
%%% EXPERIMENTAL OR NOT YET IMPLEMENTED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write.uint8_segmented = 0; % experimental: threshold segmentation for histograms with 2 distinct peaks: __/\_/\__

ADD_DEFAULT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETER / DATA SETS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2017 Maysyn22_77L_Mg5Gd_8w_a
raw_path = '/asap3/petra3/gpfs/p05/2017/data/11003950/raw/';

scan_path = [raw_path 'syn01_48L_PEEK_12w_b'];ADD
scan_path = [raw_path 'syn01_48L_PEEK_12w_c'];ADD
scan_path = [raw_path 'syn02_46R_PEEK_12w_a'];ADD
scan_path = [raw_path 'syn02_46R_PEEK_12w_b'];ADD
scan_path = [raw_path 'syn03_12R_PEEK_12w_a'];ADD
scan_path = [raw_path 'syn03_12R_PEEK_12w_b'];ADD
scan_path = [raw_path 'syn04_30R_PEEK_8w_a'];ADD
scan_path = [raw_path 'syn04_30R_PEEK_8w_b'];ADD
scan_path = [raw_path 'syn05_41R_PEEK_12w_a'];ADD
scan_path = [raw_path 'syn11_53R_Mg5Gd_12w_load_broken'];ADD

% load: att reco
raw_roi = [287 1486];
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_00'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_02'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_04'];
ref_range = [1:135, 137:162];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_06'];
ref_range = 1;ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_08'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_10'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_12'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_14'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_16'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_18'];ADD
% empty: scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_20'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_22'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_24_noload'];ADD

% load: phase reco
raw_roi = [287 1486];
do_phase_retrieval = 1;
write_8bit_segmented = 0;
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_00'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_02'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_04'];
ref_range = [1:135, 137:162];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_06'];
ref_range = 1;ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_08'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_10'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_12'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_14'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_16'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_18'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_22'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_24_noload'];ADD


raw_roi = [1211 2410];
scan_path = [raw_path 'syn14_48L_PEEK_12w_a'];ADD
scan_path = [raw_path 'syn14_48L_PEEK_12w_b'];ADD
scan_path = [raw_path 'syn15_29R_PEEK_8w_a'];ADD
scan_path = [raw_path 'syn15_29R_PEEK_8w_b'];ADD
scan_path = [raw_path 'syn16_43R_PEEK_12w_a'];ADD
scan_path = [raw_path 'syn16_43R_PEEK_12w_b'];ADD
scan_path = [raw_path 'syn17_25R_PEEK_8w_a'];ADD('r')

%% Radiography after load increase before tomography
raw_roi = [];
do_tomo = 0;
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_01'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_03'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_05'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_07'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_09'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_11'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_13'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_15'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_17'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_19'];ADD
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_21'];ADD('r')

%% Tomography CPD straw 2
raw_roi = -8;
scan_path = [raw_path 'syn22_77L_Mg5Gd_8w_a'];ADD
scan_path = [raw_path 'syn22_77L_Mg5Gd_8w_b'];ADD
scan_path = [raw_path 'syn22_80L_Mg5Gd_8w_a'];ADD
scan_path = [raw_path 'syn22_80L_Mg5Gd_8w_b'];ADD
scan_path = [raw_path 'syn22_87L_Mg5Gd_4w_a'];ADD
scan_path = [raw_path 'syn22_87L_Mg5Gd_4w_b'];ADD
scan_path = [raw_path 'syn22_88R_Mg5Gd_4w_a'];ADD
scan_path = [raw_path 'syn22_88R_Mg5Gd_4w_b'];ADD
scan_path = [raw_path 'syn22_99L_Mg5Gd_4w_a'];ADD
scan_path = [raw_path 'syn22_99L_Mg5Gd_4w_b'];ADD('r')

% CPD straw II:top
scan_path = [raw_path 'syn23_28R_PEEK_8w_a'];ADD
scan_path = [raw_path 'syn23_28R_PEEK_8w_b'];ADD('r')


%% TEST SECTION
scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_00'];
raw_roi = [287 1486];
raw_bin = 1;
do_phase_retrieval = 1;
phase_retrieval_method = 'tie';
phase_retrieval_reg_par = 2.5; 
phase_retrieval_bin_filt = 0.15; 
phase_retrieval_cutoff_frequ = 1 * pi; 
phase_padding = 1; 
write_8bit_binned = 1;
ADD

phase_retrieval_reg_par = 1.5; ADD
phase_retrieval_reg_par = 3.5; ADD
phase_retrieval_method = 'qp';phase_retrieval_reg_par = 2.5; ADD
phase_retrieval_method = 'qpcut';ADD('r')

scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_00'];
raw_roi = [287 1486];
write_float = 1; 
write_float_binned = 1; 
write_16bit = 1;
write_16bit_binned = 1;
write_8bit = 1;
write_8bit_binned = 1;
compression_method = 'histo';
compression_parameter = [0.05 0.05];
parfolder = 'test';
ADD('r')

scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_00'];
raw_roi = [287 1486];
raw_bin = 1;
butterworth_filter = 0;
write_flatcor = 0;
write_float = 1; 
write_float_binned = 1; 
write_16bit = 0;
write_16bit_binned = 0;
write_8bit = 0;
write_8bit_binned = 0;
parfolder = 'noButterworth';
ADD

scan_path = [raw_path 'syn13_55L_Mg10Gd_12w_load_00'];
raw_roi = [287 1486];
raw_bin = 1;
correlation_method =  'ssim';
write_flatcor = 0;
write_float = 1; 
write_float_binned = 1; 
write_16bit = 0;
write_16bit_binned = 0;
write_8bit = 0;
write_8bit_binned = 0;
parfolder = 'test_ssim';
ADD

correlation_method =  'ssim-ml';
parfolder = 'test_ssim-ml';
ADD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p05_reco_loop( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
