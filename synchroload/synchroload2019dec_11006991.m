function synchroload2019dec_11006991( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
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
% PARAMETER section with a 'SET_DEFAULT' statement, then the first call of
% 'ADD' will define default parameters.
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
% Created on 12-Mar-2020 by moosmanj

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

%%% SCAN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
read_flatcor = 0; % read preprocessed flatfield-corrected projections. CHECK if negative log has to be taken!
read_flatcor_path = '/asap3/petra3/gpfs/p05/2019/data/11007890/processed/AZ91_C500_ae_0/flat_corrected/rawBin2'; % subfolder of 'flat_corrected' containing projections
read_flatcor_trafo = @(im) im; %fliplr( im ); % anonymous function applied to the image which is read e.g. @(x) rot90(x)
read_sino = 0; % read preprocessed sinograms. CHECK if negative log has to be taken!
read_sino_folder = '';'test03'; % subfolder to scan path
read_sino_trafo = @(x) (x);%rot90(x); % anonymous function applied to the image which is read e.g. @(x) rot90(x)
energy = 30000; %[]; % in eV! if empty: read from log file
sample_detector_distance = [];%0.14;%[]; % in m. if empty: read from log file
eff_pixel_size = [];%0.953e-6;[]; % in m. if empty: read from log lfile. effective pixel size =  detector pixel size / magnification
pixel_scaling = 1; % to account for beam divergence if pixel size was determined (via MTF) at the wrong distance
%%% PREPROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_roi = []; % vertical and/or horizontal ROI; (1,1) coordinate = top left pixel; supports absolute, relative, negative, and mixed indexing.
% []: use full image;
% [y0 y1]: vertical ROI, skips first raw_roi(1)-1 lines, reads until raw_roi(2); if raw_roi(2) < 0 reads until end - |raw_roi(2)|; relative indexing similar.
% [y0 y1 x0 x1]: vertical + horzontal ROI, each ROI as above
% if -1: auto roi, selects vertical ROI automatically. Use only for DCM. Not working for *.raw data where images are flipped and DMM data.
% if < -1: Threshold is set as min(proj(:,:,[1 end])) + abs(raw_roi)*median(dark(:)). raw_roi=-1 defaults to min(proj(:,:,[1 end])) + 4*median(dark(:))
raw_bin = 4; % projection binning factor: integer
im_trafo = '' ;%'rot90(im,-1)'; % string to be evaluated after reading data in the case the image is flipped/rotated/etc due to changes at the beamline, e.g. 'rot90(im)'
% STITCHING/CROPPING only for scans without lateral movment. Legacy support
proj_range = []; % range of projections to be used (from all found, if empty or 1: all, if scalar: stride, if range: start:incr:end
ref_range = []; % range of flat fields to be used (from all found), if empty or 1: all. if scalar: stride, if range: start:incr:end
pixel_filter_threshold_dark = [0.01 0.005]; % Dark fields: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_threshold_flat = [0.01 0.005]; % Flat fields: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_threshold_proj = [0.01 0.005]; % Raw projection: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_radius = [3 3]; % Increase only if blobs of zeros or other artefacts are expected. Can increase processing time heavily.
ring_current_normalization = 1; % normalize flat fields and projections by ring current
image_correlation.method = 'ssim-ml';'entropy';'none';'ssim';'ssim-g';'std';'cov';'corr';'diff1-l1';'diff1-l2';'diff2-l1';'diff2-l2';'cross-entropy-12';'cross-entropy-21';'cross-entropy-x';
image_correlation.force_calc = 0; % bool. force compuation of correlation even though a (previously computed) corrlation matrix exists
image_correlation.num_flats = 3; % number of best maching flat fields used for correction
image_correlation.area_width = [1 100];%[-100 1];% correlation area: index vector or relative/absolute position of [first pix, last pix], negative indexing is supported
image_correlation.area_height = [0.2 0.8]; % correlation area: index vector or relative/absolute position of [first pix, last pix]
norm_sino = 1; % angle-wise offset-normalization of sino
strong_abs_thresh = 1; % if 1: does nothing, if < 1: flat-corrected values below threshold are set to one. Try with algebratic reco techniques.
%%% PHASE RETRIEVAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phase_retrieval.apply = 0; % See 'PhaseFilter' for detailed description of parameters !
phase_retrieval.apply_before = 0; % before stitching, interactive mode, etc. For phase-contrast data with an excentric rotation axis phase retrieval should be done afterwards. To find the rotataion axis position use this option in a first run, and then turn it of afterwards.
phase_retrieval.post_binning_factor = 1; % Binning factor after phase retrieval, but before tomographic reconstruction
phase_retrieval.method = 'qp';'tie';'qpcut';'dpc'; %'qp' 'ctf' 'tie' 'qp2' 'qpcut'
phase_retrieval.reg_par = 2.5; % regularization parameter. larger values tend to blurrier images. smaller values tend to original data.
phase_retrieval.bin_filt = 0.1; % threshold for quasiparticle retrieval 'qp', 'qp2'
phase_retrieval.cutoff_frequ = 2 * pi; % in radian. frequency cutoff in Fourier space for 'qpcut' phase retrieval
phase_retrieval.padding = 1; % padding of intensities before phase retrieval, 0: no padding
dpc_steps = 5;
dpc_bin = 8;
%%% TOMOGRAPHY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tomo.run = 0; % run tomographic reconstruction
tomo.run_interactive_mode = 1; % if tomo.run = 0, use to determine rot axis positions without processing the full tomogram;
tomo.reco_mode = '3D';'slice';  % slice-wise or full 3D backprojection. 'slice': volume must be centered at origin & no support of rotation axis tilt, reco binning, save compressed
tomo.vol_size = []; %[-.5 .5 -.5 .5 -0.5 0.5];% 6-component vector [xmin xmax ymin ymax zmin zmax], for excentric rot axis pos / extended FoV;. if empty, volume is centerd within tomo.vol_shape. unit voxel size is assumed. if smaller than 10 values are interpreted as relative size w.r.t. the detector size. Take care bout minus signs! Note that if empty vol_size is dependent on the rotation axis position.
tomo.vol_shape = []; %[1 1 1] shape (# voxels) of reconstruction volume. used for excentric rot axis pos. if empty, inferred from 'tomo.vol_size'. in absolute numbers of voxels or in relative number w.r.t. the default volume which is given by the detector width and height.
tomo.rot_angle_full_range = [];% 2 * pi ;[]; % in radians. if []: full angle of rotation including additional increment, or array of angles. if empty full rotation angles is determined automatically to pi or 2 pi
tomo.rot_angle_offset = pi; % global rotation of reconstructed volume
tomo.rot_axis_offset = 1.4;% -230.2;[]; % rotation axis offset w.r.t to the image center. Assuming the rotation axis position to be centered in the FOV for standard scan, the offset should be close to zero.
tomo.rot_axis_position = []; % if empty use automatic computation. EITHER OFFSET OR POSITION MUST BE EMPTY. YOU MUST NOT USE BOTH!
tomo.rot_axis_offset_shift = []; %[]; % absolute lateral movement in pixels during fly-shift-scan, overwrite lateral shift read out from hdf5 log
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
tomo.algorithm = 'fbp';'sirt'; 'cgls';'sart';'em';'fbp-astra'; % SART/EM only work for 3D reco mode
tomo.iterations = 40; % for 'sirt' or 'cgls'.
tomo.MinConstraint = []; % If specified, all values below MinConstraint will be set to MinConstraint. This can be used to enforce non-negative reconstructions, for example.
tomo.MaxConstraint = []; % If specified, all values above MaxConstraint will be set to MaxConstraint.
%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write.path = ''; %'/gpfs/petra3/scratch/moosmanj'; % absolute path were output data will be stored. !!overwrites the write.to_scratch flag. if empty uses the beamtime directory and either 'processed' or 'scratch_cc'
write.to_scratch = 0; % write to 'scratch_cc' instead of 'processed'
write.deleteFiles = 0; % delete files already existing in output folders. Useful if number or names of files differ when reprocessing.
write.beamtimeID = ''; % string (regexp),typically beamtime ID, mandatory if 'write.deleteFiles' is true (safety check)
write.scan_name_appendix = ''; % appendix to the output folder name which defaults to the scan name
write.parfolder = '';% parent folder to 'reco', 'sino', 'phase', and 'flat_corrected'
write.subfolder_flatcor = ''; % subfolder in 'flat_corrected'
write.subfolder_phase_map = ''; % subfolder in 'phase_map'
write.subfolder_sino = ''; % subfolder in 'sino'
write.subfolder_reco = ''; % subfolder in 'reco'
write.flatcor = 1; % save preprocessed flat corrected projections
write.phase_map = 1; % save phase maps (if phase retrieval is not 0)
write.sino = 0; % save sinograms (after preprocessing & before FBP filtering and phase retrieval)
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
write.uint8_segmented = 0; % experimental: threshold segmentaion for histograms with 2 distinct peaks: __/\_/\__
%%% INTERACTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.visual_output = 0; % show images and plots during reconstruction
interactive_mode.rot_axis_pos = 0; % reconstruct slices with dif+ferent rotation axis offsets
interactive_mode.rot_axis_pos_default_search_range = []; % if empty: asks for search range when entering interactive mode
interactive_mode.rot_axis_tilt = 0; % reconstruct slices with different offset AND tilts of the rotation axis
interactive_mode.rot_axis_tilt_default_search_range = []; % if empty: asks for search range when entering interactive mode
interactive_mode.lamino = 0; % find laminography tilt instead camera tilt
interactive_mode.angles = 0; % reconstruct slices with different scalings of angles
interactive_mode.angle_scaling_default_search_range = []; % if empty: use a variaton of -/+5 * (angle increment / maximum angle)
interactive_mode.slice_number = 0.5; % default slice number. if in [0,1): relative, if in (1, N]: absolute
interactive_mode.phase_retrieval = 0; % Interactive retrieval to determine regularization parameter
interactive_mode.phase_retrieval_default_search_range = []; % if empty: asks for search range when entering interactive mode, otherwise directly start with given search range
%%% HARDWARE / SOFTWARE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.use_cluster = 0; % if available: on MAXWELL nodes disp/nova/wga/wgs cluster computation can be used. Recommended only for large data sets since parpool creation and data transfer implies a lot of overhead.
par.poolsize = 0.75; % number of workers used in a local parallel pool. if 0: use current config. if >= 1: absolute number. if 0 < poolsize < 1: relative amount of all cores to be used. if SLURM scheduling is available, a default number of workers is used.
par.poolsize_gpu_limit_factor = 0.7; % Relative amount of GPU memory used for preprocessing during parloop. High values speed up Proprocessing, but increases out-of-memory failure
tomo.astra_link_data = 1; % ASTRA data objects become references to Matlab arrays. Reduces memory issues.
tomo.astra_gpu_index = [2 4:6];
par.gpu_index = tomo.astra_gpu_index;
par.use_gpu_in_parfor = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set default. Defines default parameter set to be used.
SET_DEFAULT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETER / DATA SETS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_path = '/asap3/petra3/gpfs/p07/2019/data/11006991/raw/';

scan_path = [raw_path 'hzg_ind_01_cork_a']; ADD

scan_path = [raw_path 'hzg_ind_02_silicon_a']; ADD

scan_path = [raw_path 'hzg_ind_02_silicon_b']; ADD

scan_path = [raw_path 'hzg_p07_2019_11']; ADD

scan_path = [raw_path 'syn001_Ti_12w_47R']; ADD

scan_path = [raw_path 'syn002_Ti_12w_47R_z0140_t130ms']; ADD
scan_path = [raw_path 'syn004_Ti_12w_47R_z0280_t130ms']; ADD

%% syn005 1ms
scan_path = [raw_path 'syn005_Ti_12w_47R_z0030_t001ms']; 
tomo.rot_axis_offset = 0.4 * 4 / raw_bin;
tomo.reco_mode = 'slice';

tomo.algorithm = 'fbp';
write.subfolder_reco = 'fbp';
ADD

tomo.algorithm = 'sirt';
tomo.iterations = 200;
write.subfolder_reco = sprintf( '%s%u',  tomo.algorithm, tomo.iterations );
ADD

tomo.algorithm = 'cgls';
write.subfolder_reco = sprintf( '%s%u',  tomo.algorithm, tomo.iterations );
ADD

tomo.algorithm = 'sart';
write.subfolder_reco = sprintf( '%s%u',  tomo.algorithm, tomo.iterations );
ADD

tomo.algorithm = 'em';
write.subfolder_reco = sprintf( '%s%u',  tomo.algorithm, tomo.iterations );
ADD

%% syn006 5ms
scan_path = [raw_path 'syn006_Ti_12w_47R_z0030_t005ms'];

tomo.algorithm = 'fbp';
write.subfolder_reco = sprintf( '%s',  tomo.algorithm  );
ADD

tomo.iterations = 200;

tomo.algorithm = 'sirt';
write.subfolder_reco = sprintf( '%s%u',  tomo.algorithm, tomo.iterations );
ADD

tomo.algorithm = 'cgls';
write.subfolder_reco = sprintf( '%s%u',  tomo.algorithm, tomo.iterations );
ADD

tomo.algorithm = 'sart';
write.subfolder_reco = sprintf( '%s%u',  tomo.algorithm, tomo.iterations );
ADD

tomo.algorithm = 'em';
write.subfolder_reco = sprintf( '%s%u',  tomo.algorithm, tomo.iterations );
ADD

tomo.astra_gpu_index = [];
par.gpu_index = tomo.astra_gpu_index;

%% syn007 10ms
tomo.algorithm = 'fbp';
write.subfolder_reco = sprintf( '%s%u',  tomo.algorithm );
scan_path = [raw_path 'syn007_Ti_12w_47R_z0030_t010ms']; ADD

%% syn008 135ms
scan_path = [raw_path 'syn008_Ti_12w_47R_z0030_t135ms'];

tomo.reco_mode = 'slice';

tomo.algorithm = 'fbp';
write.subfolder_reco = sprintf( '%s%u',  tomo.algorithm );
ADD

tomo.algorithm = 'sirt';
write.subfolder_reco = sprintf( '%s%u',  tomo.algorithm, tomo.iterations );
ADD

tomo.algorithm = 'cgls';
write.subfolder_reco = sprintf( '%s%u',  tomo.algorithm, tomo.iterations );
ADD

%%
tomo.algorithm = 'fbp';
write.subfolder_reco = 'fbp';
tomo.rot_axis_offset = 0.0 * 4 / raw_bin;
scan_path = [raw_path 'syn009_Ti_12w_47R_z0140_t001ms']; ADD
scan_path = [raw_path 'syn010_Ti_12w_47R_z0140_t005ms']; ADD
scan_path = [raw_path 'syn011_Ti_12w_47R_z0140_t010ms']; ADD
scan_path = [raw_path 'syn012_Ti_12w_47R_z0140_t135ms']; ADD

tomo.rot_axis_offset = -1.7 * 4 / raw_bin;
scan_path = [raw_path 'syn013_Ti_12w_47R_z0280_t001ms']; ADD
scan_path = [raw_path 'syn014_Ti_12w_47R_z0280_t005ms']; ADD
scan_path = [raw_path 'syn015_Ti_12w_47R_z0280_t010ms']; ADD
scan_path = [raw_path 'syn016_Ti_12w_47R_z0280_t135ms']; ADD

tomo.rot_axis_offset = -4.6 * 4 / raw_bin;
scan_path = [raw_path 'syn017_Ti_12w_47R_z0560_t001ms']; ADD
scan_path = [raw_path 'syn018_Ti_12w_47R_z0560_t005ms']; ADD
scan_path = [raw_path 'syn019_Ti_12w_47R_z0560_t010ms']; ADD
scan_path = [raw_path 'syn020_Ti_12w_47R_z0560_t135ms']; ADD

tomo.rot_axis_offset = -11 * 4 / raw_bin;
scan_path = [raw_path 'syn021_Ti_12w_47R_z1120_t001ms']; ADD
scan_path = [raw_path 'syn022_Ti_12w_47R_z1120_t005ms']; ADD
scan_path = [raw_path 'syn023_Ti_12w_47R_z1120_t010ms']; ADD
scan_path = [raw_path 'syn024_Ti_12w_47R_z1120_t135ms']; ADD

tomo.rot_axis_offset = -8.35 * 4 / raw_bin;
scan_path = [raw_path 'syn025_Ti_12w_47R_z0900_t001ms_dpc']; ADD
scan_path = [raw_path 'syn026_Ti_12w_47R_z0900_t005ms_dpc']; ADD
scan_path = [raw_path 'syn027_Ti_12w_47R_z0900_t010ms_dpc']; ADD

scan_path = [raw_path 'syn31_621732li_6w_Ti_prox1_01']; ADD
scan_path = [raw_path 'syn31_621732li_6w_Ti_prox1_02']; ADD
scan_path = [raw_path 'syn31_621732li_6w_Ti_prox1_03']; ADD
scan_path = [raw_path 'syn31_621732li_6w_Ti_prox1_04']; ADD
scan_path = [raw_path 'syn31_621732li_6w_Ti_prox1_05']; ADD
scan_path = [raw_path 'syn31_621732li_6w_Ti_prox1_06']; ADD
scan_path = [raw_path 'syn31_621732li_6w_Ti_prox1_07']; ADD
scan_path = [raw_path 'syn31_621732li_6w_Ti_prox1_08']; ADD
scan_path = [raw_path 'syn31_621732li_6w_Ti_prox1_09']; ADD
scan_path = [raw_path 'syn31_621732li_6w_Ti_prox1_10']; ADD
scan_path = [raw_path 'syn31_621732li_6w_Ti_prox1_11']; ADD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p05_reco_loop( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
