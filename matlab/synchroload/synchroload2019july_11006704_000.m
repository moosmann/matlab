function p05_reco_loop_synchroload2019july_11006704_000( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
% Created on 05-Aug-2019 by moosmanj

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

fast_reco = 1; % !!! OVERWRITES SOME PARAMETERS SET BELOW !!!
stop_after_data_reading(1) = 0; % for data analysis, before flat field correlation
stop_after_proj_flat_correlation(1) = 0; % for data analysis, after flat field correlation

%%% SCAN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scan_path = pwd; % sting/pwd. pwd: change to directory of the scan to be reconstructed, sting: absolute scan path
read_flatcor = 0; % read preprocessed flatfield-corrected projections. CHECK if negative log has to be taken!
read_flatcor_path = ''; % subfolder of 'flat_corrected' containing projections
read_sino = 0; % read preprocessed sinograms. CHECK if negative log has to be taken!
read_sino_folder = ''; % subfolder to scan path
energy = []; % in eV! if empty: read from log file
sample_detector_distance = []; % in m. if empty: read from log file
eff_pixel_size = []; % in m. if empty: read from log file. effective pixel size =  detector pixel size / magnification
pix_scaling = 1; % to account for beam divergence if pixel size was determined (via MTF) at the wrong distance
%%% PREPROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_roi = []; % vertical and/or horizontal ROI; (1,1) coordinate = top left pixel; supports absolute, relative, negative, and mixed indexing.
% []: use full image;
% [y0 y1]: vertical ROI, skips first raw_roi(1)-1 lines, reads until raw_roi(2); if raw_roi(2) < 0 reads until end - |raw_roi(2)|; relative indexing similar.
% [y0 y1 x0 x1]: vertical + horzontal ROI, each ROI as above
% if -1: auto roi, selects vertical ROI automatically. Use only for DCM. Not working for *.raw data where images are flipped and DMM data.
% if < -1: Threshold is set as min(proj(:,:,[1 end])) + abs(raw_roi)*median(dark(:)). raw_roi=-1 defaults to min(proj(:,:,[1 end])) + 4*median(dark(:))
raw_bin = 4; % projection binning factor: integer
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
pixel_filter_radius = [3 3]; % Increase only if blobs of zeros or other artefacts are expected. Can increase processing time heavily.
ring_current_normalization = 1; % normalize flat fields and projections by ring current
image_correlation.method = 'ssim-ml';'entropy';'diff';'shift';'ssim';'std';'cov';'corr';'cross-entropy12';'cross-entropy21';'cross-entropyx';
image_correlation.force_calc = 0; % bool. force compuation of correlation even though a (previously computed) corrlation matrix exists
image_correlation.num_flats = 1; % number of flat fields used for average/median of flats. for 'shift'-correlation its the maximum number
image_correlation.area_width = [0 0.01];%[0.98 1];% correlation area: index vector or relative/absolute position of [first pix, last pix]
image_correlation.area_height = [0.3 0.7]; % correlation area: index vector or relative/absolute position of [first pix, last pix]
ring_filter.apply = 0; % ring artifact filter (use only for scans without lateral sample movement)
ring_filter.apply_before_stitching = 0; % ! Consider when phase retrieval is applied !
ring_filter.method = 'jm'; 'wavelet-fft';
ring_filter.waveletfft.dec_levels = 2:5; % decomposition levels for 'wavelet-fft'
ring_filter.waveletfft.wname = 'db25';'db30'; % wavelet type, see 'FilterStripesCombinedWaveletFFT' or 'waveinfo'
ring_filter.waveletfft.sigma = 2.4; %  suppression factor for 'wavelet-fft'
ring_filter.jm.median_width = 11; % multiple widths are applied consecutively, eg [3 11 21 31 39];
strong_abs_thresh = 1; % if 1: does nothing, if < 1: flat-corrected values below threshold are set to one. Try with algebratic reco techniques.
%%% PHASE RETRIEVAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phase_retrieval.apply = 0; % See 'PhaseFilter' for detailed description of parameters !
phase_retrieval.apply_before = 0; % before stitching, interactive mode, etc. For phase-contrast data with an excentric rotation axis phase retrieval should be done afterwards. To find the rotataion axis position use this option in a first run, and then turn it of afterwards.
phase_retrieval.post_binning_factor = 1; % Binning factor after phase retrieval, but before tomographic reconstruction
phase_retrieval.method = 'tie';'qp';'qpcut'; %'qp' 'ctf' 'tie' 'qp2' 'qpcut'
phase_retrieval.reg_par = 2; % regularization parameter. larger values tend to blurrier images. smaller values tend to original data.
phase_retrieval.bin_filt = 0.1; % threshold for quasiparticle retrieval 'qp', 'qp2'
phase_retrieval.cutoff_frequ = 2 * pi; % in radian. frequency cutoff in Fourier space for 'qpcut' phase retrieval
phase_retrieval.padding = 1; % padding of intensities before phase retrieval, 0: no padding
%%% TOMOGRAPHY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tomo.run = 1; % run tomographic reconstruction
tomo.run_interactive_mode = 1; % if tomo.run = 0, intendet to determine rot axis positions without processing the full tomogram;
tomo.reco_mode = '3D'; 'slice'; % slice-wise or full 3D backprojection. 'slice': volume must be centered at origin & no support of rotation axis tilt, reco binning, save compressed
tomo.vol_size = [];%[-1 1 -1 1 -0.5 0.5];% 6-component vector [xmin xmax ymin ymax zmin zmax], for excentric rot axis pos / extended FoV;. if empty, volume is centerd within tomo.vol_shape. unit voxel size is assumed. if smaller than 10 values are interpreted as relative size w.r.t. the detector size. Take care bout minus signs! Note that if empty vol_size is dependent on the rotation axis position.
tomo.vol_shape = []; %[1 1 1] shape (# voxels) of reconstruction volume. used for excentric rot axis pos. if empty, inferred from 'tomo.vol_size'. in absolute numbers of voxels or in relative number w.r.t. the default volume which is given by the detector width and height.
tomo.rot_angle.full_range = []; % in radians. if []: full angle of rotation including additional increment, or array of angles. if empty full rotation angles is determined automatically to pi or 2 pi
tomo.rot_angle.offset = pi; % global rotation of reconstructed volume
tomo.rot_axis.offset = 0; % if empty use automatic computation
tomo.rot_axis.position = []; % if empty use automatic computation. EITHER OFFSET OR POSITION MUST BE EMPTY. YOU MUST NOT USE BOTH!
tomo.rot_axis.offset_shift_range = []; %[]; % absolute lateral movement in pixels during fly-shift-scan, overwrite lateral shift read out from hdf5 log
tomo.rot_axis.tilt = 0; % in rad. camera tilt w.r.t rotation axis. if empty calculate from registration of projections at 0 and pi
tomo.rot_axis.corr_area1 = []; % ROI to correlate projections at angles 0 & pi. Use [0.75 1] or so for scans with an excentric rotation axis 
tomo.rot_axis.corr_area2 = []; % ROI to correlate projections at angles 0 & pi 
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
write.deleteFiles = 0; % delete files already existing in output folders. Useful if number or names of files differ when reprocessing.
write.beamtimeID = ''; % string (regexp),typically beamtime ID, mandatory if 'write.deleteFiles' is true (safety check)
write.scan_name_appendix = ''; % appendix to the output folder name which defaults to the scan name
write.parfolder = '';% parent folder to 'reco', 'sino', 'phase', and 'flat_corrected'
write.subfolder.flatcor = ''; % subfolder in 'flat_corrected'
write.subfolder.phase_map = ''; % subfolder in 'phase_map'
write.subfolder.sino = ''; % subfolder in 'sino'
write.subfolder.reco = ''; % subfolder in 'reco'
write.flatcor = 0; % save preprocessed flat corrected projections
write.flatcor_shift_cropped = 0; % save lateral shift corrected projections, projections are not interpolated, but cropped to nearest integer pixel
write.phase_map = 0; % save phase maps (if phase retrieval is not 0)
write.sino = 0; % save sinograms (after preprocessing & before FBP filtering and phase retrieval)
write.sino_shift_cropped = 0; % save cropped sinos without lateral shift
write.phase_sino = 0; % save sinograms of phase maps
write.reco = 1; % save reconstructed slices (if tomo.run=1)
write.float = 1; % single precision (32-bit float) tiff
write.uint16 = 0; % save 16bit unsigned integer tiff using 'write.compression.method'
write.uint8 = 0; % save binned 8bit unsigned integer tiff using 'write.compression.method'
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
interactive_mode.rot_axis_pos_default_search_range = []; % if empty: asks for search range when entering interactive mode, otherwise directly start with given search range
interactive_mode.rot_axis_tilt = 0; % reconstruct slices with different offset AND tilts of the rotation axis
interactive_mode.lamino = 0; % find laminography tilt instead camera rotation
interactive_mode.fixed_other_tilt = 0; % fixed other tilt
interactive_mode.slice_number = 0.5; % default slice number. if in [0,1): relative, if in (1, N]: absolute
interactive_mode.phase_retrieval = 1; % Interactive retrieval to determine regularization parameter
interactive_mode.phase_retrieval_default_search_range = []; % if empty: asks for search range when entering interactive mode, otherwise directly start with given search range
%%% HARDWARE / SOFTWARE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
use_cluster = 0; % if available: on MAXWELL nodes disp/nova/wga/wgs cluster computation can be used. Recommended only for large data sets since parpool creation and data transfer implies a lot of overhead.
poolsize = 0.8; % number of workers used in a local parallel pool. if 0: use current config. if >= 1: absolute number. if 0 < poolsize < 1: relative amount of all cores to be used. if SLURM scheduling is available, a default number of workers is used.
tomo.astra_link_data = 1; % ASTRA data objects become references to Matlab arrays. Reduces memory issues.
tomo.astra_gpu_index = []; % GPU Device index to use, Matlab notation: index starts from 1. default: [], uses all
%%% EXPERIMENTAL OR NOT YET IMPLEMENTED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write.uint8_segmented = 0; % experimental: threshold segmentation for histograms with 2 distinct peaks: __/\_/\__

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set default. Defines default parameter set to be used.
SET_DEFAULT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETER / DATA SETS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_path = '/asap3/petra3/gpfs/p05/2019/data/11006704/raw/';

scan_path = [raw_path 'syn001_35R_Ti_8w_000']; ADD
scan_path = [raw_path 'syn001_35R_Ti_8w_001']; ADD
scan_path = [raw_path 'syn001_35R_Ti_8w_002']; ADD
scan_path = [raw_path 'syn001_35R_Ti_8w_003']; ADD
scan_path = [raw_path 'syn001_35R_Ti_8w_004']; ADD
scan_path = [raw_path 'syn001_35R_Ti_8w_005']; ADD
scan_path = [raw_path 'syn001_35R_Ti_8w_006']; ADD
scan_path = [raw_path 'syn001_35R_Ti_8w_007']; ADD
scan_path = [raw_path 'syn001_35R_Ti_8w_008']; ADD
scan_path = [raw_path 'syn001_35R_Ti_8w_009']; ADD

scan_path = [raw_path 'syn002_65L_Mg5Gd_12w_000']; ADD
scan_path = [raw_path 'syn002_65L_Mg5Gd_12w_001']; ADD
scan_path = [raw_path 'syn002_65L_Mg5Gd_12w_002']; ADD
scan_path = [raw_path 'syn002_65L_Mg5Gd_12w_003']; ADD
scan_path = [raw_path 'syn002_65L_Mg5Gd_12w_004']; ADD
scan_path = [raw_path 'syn002_65L_Mg5Gd_12w_005']; ADD
scan_path = [raw_path 'syn002_65L_Mg5Gd_12w_006']; ADD
scan_path = [raw_path 'syn002_65L_Mg5Gd_12w_007']; ADD
scan_path = [raw_path 'syn002_65L_Mg5Gd_12w_008']; ADD
scan_path = [raw_path 'syn002_65L_Mg5Gd_12w_009']; ADD

scan_path = [raw_path 'syn003_47R_Ti_12w_000']; ADD
scan_path = [raw_path 'syn003_47R_Ti_12w_001']; ADD
scan_path = [raw_path 'syn003_47R_Ti_12w_002']; ADD
scan_path = [raw_path 'syn003_47R_Ti_12w_003']; ADD
scan_path = [raw_path 'syn003_47R_Ti_12w_004']; ADD
scan_path = [raw_path 'syn003_47R_Ti_12w_005']; ADD
scan_path = [raw_path 'syn003_47R_Ti_12w_006']; ADD
scan_path = [raw_path 'syn003_47R_Ti_12w_008']; ADD
scan_path = [raw_path 'syn003_47R_Ti_12w_012']; ADD
scan_path = [raw_path 'syn003_47R_Ti_12w_013']; ADD

scan_path = [raw_path 'syn004_71L_Mg5Gd_12w_000']; ADD
scan_path = [raw_path 'syn004_71L_Mg5Gd_12w_001']; ADD
scan_path = [raw_path 'syn004_71L_Mg5Gd_12w_002']; ADD
scan_path = [raw_path 'syn004_71L_Mg5Gd_12w_003']; ADD
scan_path = [raw_path 'syn004_71L_Mg5Gd_12w_004']; ADD
scan_path = [raw_path 'syn004_71L_Mg5Gd_12w_005']; ADD
scan_path = [raw_path 'syn004_71L_Mg5Gd_12w_006']; ADD
scan_path = [raw_path 'syn004_71L_Mg5Gd_12w_007']; ADD
scan_path = [raw_path 'syn004_71L_Mg5Gd_12w_008']; ADD
scan_path = [raw_path 'syn004_71L_Mg5Gd_12w_009']; ADD
scan_path = [raw_path 'syn004_71L_Mg5Gd_12w_010']; ADD
scan_path = [raw_path 'syn004_71L_Mg5Gd_12w_011']; ADD
scan_path = [raw_path 'syn004_71L_Mg5Gd_12w_011_old']; ADD
scan_path = [raw_path 'syn004_71L_Mg5Gd_12w_012']; ADD
scan_path = [raw_path 'syn004_71L_Mg5Gd_12w_013']; ADD
scan_path = [raw_path 'syn004_71L_Mg5Gd_12w_014']; ADD
scan_path = [raw_path 'syn004_71L_Mg5Gd_12w_015']; ADD

scan_path = [raw_path 'syn005_40L_Peek_12w_000']; ADD
scan_path = [raw_path 'syn007_47L_Peek_12w_000']; ADD
scan_path = [raw_path 'syn007_47L_Peek_12w_001']; ADD
scan_path = [raw_path 'syn007_47L_Peek_12w_002']; ADD
scan_path = [raw_path 'syn007_47L_Peek_12w_003']; ADD
scan_path = [raw_path 'syn007_47L_Peek_12w_004']; ADD
scan_path = [raw_path 'syn007_47L_Peek_12w_005']; ADD
scan_path = [raw_path 'syn007_47L_Peek_12w_006']; ADD
scan_path = [raw_path 'syn007_47L_Peek_12w_007']; ADD
scan_path = [raw_path 'syn007_47L_Peek_12w_008']; ADD
scan_path = [raw_path 'syn007_47L_Peek_12w_009']; ADD
scan_path = [raw_path 'syn007_47L_Peek_12w_010']; ADD
scan_path = [raw_path 'syn007_47L_Peek_12w_011']; ADD
scan_path = [raw_path 'syn007_47L_Peek_12w_012']; ADD
scan_path = [raw_path 'syn007_47L_Peek_12w_013']; ADD

scan_path = [raw_path 'syn008_47R_Ti_12w_001']; ADD
scan_path = [raw_path 'syn008_47R_Ti_12w_002']; ADD
scan_path = [raw_path 'syn008_47R_Ti_12w_003']; ADD
scan_path = [raw_path 'syn008_47R_Ti_12w_004']; ADD
scan_path = [raw_path 'syn008_47R_Ti_12w_005']; ADD
scan_path = [raw_path 'syn008_47R_Ti_12w_006']; ADD
scan_path = [raw_path 'syn008_47R_Ti_12w_007']; ADD
scan_path = [raw_path 'syn008_47R_Ti_12w_008']; ADD
scan_path = [raw_path 'syn008_47R_Ti_12w_009']; ADD
scan_path = [raw_path 'syn008_47R_Ti_12w_010']; ADD
scan_path = [raw_path 'syn008_47R_Ti_12w_011']; ADD
scan_path = [raw_path 'syn008_47R_Ti_12w_012']; ADD
scan_path = [raw_path 'syn008_47R_Ti_12w_013']; ADD
scan_path = [raw_path 'syn008_47R_Ti_12w_014']; ADD
scan_path = [raw_path 'syn008_47R_Ti_12w_015']; ADD
scan_path = [raw_path 'syn008_47R_Ti_12w_016']; ADD
scan_path = [raw_path 'syn008_47R_Ti_12w_017']; ADD
scan_path = [raw_path 'syn008_47R_Ti_12w_018']; ADD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p05_reco_loop( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
