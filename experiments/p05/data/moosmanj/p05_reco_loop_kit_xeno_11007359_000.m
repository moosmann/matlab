function p05_reco_loop_kit_xeno_11007359_000( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
% Created on 20-Dec-2019 by moosmanj

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

fast_reco.run = 0; 

%%% SCAN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
read_flatcor = 0; % read preprocessed flatfield-corrected projections. CHECK if negative log has to be taken!
read_flatcor_path = ''; % subfolder of 'flat_corrected' containing projections
read_flatcor_trafo = @(im) fliplr( im ); % anonymous function applied to the image which is read e.g. @(x) rot90(x)
read_sino = 0; % read preprocessed sinograms. CHECK if negative log has to be taken!
read_sino_folder = 'trans03'; % subfolder to scan path
read_sino_trafo = @(x) (x);%rot90(x); % anonymous function applied to the image which is read e.g. @(x) rot90(x)
energy = []; % in eV! if empty: read from log file
sample_detector_distance = []; % in m. if empty: read from log file
eff_pixel_size = []; % in m. if empty: read from log lfile. effective pixel size =  detector pixel size / magnification
pixel_scaling = 1; % to account for beam divergence if pixel size was determined (via MTF) at the wrong distance
%%% PREPROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_roi = []; % vertical and/or horizontal ROI; (1,1) coordinate = top left pixel; supports absolute, relative, negative, and mixed indexing.
% []: use full image;
% [y0 y1]: vertical ROI, skips first raw_roi(1)-1 lines, reads until raw_roi(2); if raw_roi(2) < 0 reads until end - |raw_roi(2)|; relative indexing similar.
% [y0 y1 x0 x1]: vertical + horzontal ROI, each ROI as above
% if -1: auto roi, selects vertical ROI automatically. Use only for DCM. Not working for *.raw data where images are flipped and DMM data.
% if < -1: Threshold is set as min(proj(:,:,[1 end])) + abs(raw_roi)*median(dark(:)). raw_roi=-1 defaults to min(proj(:,:,[1 end])) + 4*median(dark(:))
raw_bin = 2; % projection binning factor: integer
im_trafo = '' ;%'rot90(im,-1)'; % string to be evaluated after reading data in the case the image is flipped/rotated/etc due to changes at the beamline, e.g. 'rot90(im)'
proj_range = []; % range of projections to be used (from all found, if empty or 1: all, if scalar: stride, if range: start:incr:end
ref_range = []; % range of flat fields to be used (from all found), if empty or 1: all. if scalar: stride, if range: start:incr:end
pixel_filter_threshold_dark = [0.01 0.005]; % Dark fields: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_threshold_flat = [0.01 0.005]; % Flat fields: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_threshold_proj = [0.01 0.005]; % Raw projection: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_radius = [3 3]; % Increase only if blobs of zeros or other artefacts are expected. Can increase processing time heavily.
ring_current_normalization = 0; % normalize flat fields and projections by ring current
image_correlation.method = 'ssim-ml';'entropy';'ssim';'ssim-g';'std';'cov';'corr';'diff1-l1';'diff1-l2';'diff2-l1';'diff2-l2';'cross-entropy-12';'cross-entropy-21';'cross-entropy-x';
image_correlation.force_calc = 0; % bool. force compuation of correlation even though a (previously computed) corrlation matrix exists
image_correlation.num_flats = 3; % number of best maching flat fields used for correction
image_correlation.area_width = [1 100];%[-100 1];% correlation area: index vector or relative/absolute position of [first pix, last pix], negative indexing is supported
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
phase_retrieval.apply = 1; % See 'PhaseFilter' for detailed description of parameters !
phase_retrieval.apply_before = 0; % before stitching, interactive mode, etc. For phase-contrast data with an excentric rotation axis phase retrieval should be done afterwards. To find the rotataion axis position use this option in a first run, and then turn it of afterwards.
phase_retrieval.post_binning_factor = 1; % Binning factor after phase retrieval, but before tomographic reconstruction
phase_retrieval.method = 'tie';'qp';'qpcut'; %'qp' 'ctf' 'tie' 'qp2' 'qpcut'
phase_retrieval.reg_par = 1.7; % regularization parameter. larger values tend to blurrier images. smaller values tend to original data.
phase_retrieval.bin_filt = 0.1; % threshold for quasiparticle retrieval 'qp', 'qp2'
phase_retrieval.cutoff_frequ = 2 * pi; % in radian. frequency cutoff in Fourier space for 'qpcut' phase retrieval
phase_retrieval.padding = 1; % padding of intensities before phase retrieval, 0: no padding
%%% TOMOGRAPHY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tomo.run = 1; % run tomographic reconstruction
tomo.run_interactive_mode = 1; % if tomo.run = 0, use to determine rot axis positions without processing the full tomogram;
tomo.reco_mode = '3D'; 'slice'; % slice-wise or full 3D backprojection. 'slice': volume must be centered at origin & no support of rotation axis tilt, reco binning, save compressed
tomo.vol_size = []; %[-.5 .5 -.5 .5 -0.5 0.5];% 6-component vector [xmin xmax ymin ymax zmin zmax], for excentric rot axis pos / extended FoV;. if empty, volume is centerd within tomo.vol_shape. unit voxel size is assumed. if smaller than 10 values are interpreted as relative size w.r.t. the detector size. Take care bout minus signs! Note that if empty vol_size is dependent on the rotation axis position.
tomo.vol_shape = []; %[1 1 1] shape (# voxels) of reconstruction volume. used for excentric rot axis pos. if empty, inferred from 'tomo.vol_size'. in absolute numbers of voxels or in relative number w.r.t. the default volume which is given by the detector width and height.
tomo.rot_angle.full_range = 2*pi; []; % in radians. if []: full angle of rotation including additional increment, or array of angles. if empty full rotation angles is determined automatically to pi or 2 pi
tomo.rot_angle.offset = pi; % global rotation of reconstructed volume
tomo.rot_axis.offset = []; % rotation axis offset w.r.t to the image center. Assuming the rotation axis position to be centered in the FOV for standard scan, the offset should be close to zero.
tomo.rot_axis.position = []; % if empty use automatic computation. EITHER OFFSET OR POSITION MUST BE EMPTY. YOU MUST NOT USE BOTH!
tomo.rot_axis.offset_shift = []; %[]; % absolute lateral movement in pixels during fly-shift-scan, overwrite lateral shift read out from hdf5 log
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
write.flatcor = 1; % save preprocessed flat corrected projections
write.phase_map = 0; % save phase maps (if phase retrieval is not 0)
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
write.compression.method = 'outlier';'histo';'full'; 'std'; 'threshold'; % method to compression dynamic range into [0, 1]
write.compression.parameter = [0.02 0.02]; % compression-method specific parameter
% dynamic range is compressed s.t. new dynamic range assumes
% 'outlier' : [LOW, HIGH] = write.compression.parameter, eg. [0.01 0.03], outlier given in percent, if scalear LOW = HIGH.
% 'full' : full dynamic range is used
% 'threshold' : [LOW HIGH] = write.compression.parameter, eg. [-0.01 1]
% 'std' : NUM = write.compression.parameter, mean +/- NUM*std, dynamic range is rescaled to within -/+ NUM standard deviations around the mean value
% 'histo' : [LOW HIGH] = write.compression.parameter (100*LOW)% and (100*HIGH)% of the original histogram, e.g. [0.02 0.02]
%%% INTERACTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.visual_output = 0; % show images and plots during reconstruction
interactive_mode.rot_axis_pos = 0; % reconstruct slices with dif+ferent rotation axis offsets
interactive_mode.rot_axis_pos_default_search_range = []; % if empty: asks for search range when entering interactive mode
interactive_mode.rot_axis_tilt = 0; % reconstruct slices with different offset AND tilts of the rotation axis
interactive_mode.rot_axis_tilt_default_search_range = []; % if empty: asks for search range when entering interactive mode
interactive_mode.lamino = 0; % find laminography tilt instead camera rotation
interactive_mode.fixed_other_tilt = 0; % fixed other tilt
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
tomo.astra_gpu_index = []; % GPU Device index to use, Matlab notation: index starts from 1. default: [], uses all
%%% EXPERIMENTAL OR NOT YET IMPLEMENTED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write.uint8_segmented = 0; % experimental: threshold segmentaion for histograms with 2 distinct peaks: __/\_/\__

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set default. Defines default parameter set to be used.
SET_DEFAULT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETER / DATA SETS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_path = '/asap3/petra3/gpfs/p05/2019/data/11007359/raw/';

pixel_filter_radius = [5 5]; 

% poor image quality: scintillator extremely dirty and strong beam
% variations, almost impossible to find a proper reference image

%% Test Section

phase_retrieval.apply = 0;
par.visual_output = 1; 
write.flatcor = 1; 
interactive_mode.rot_axis_pos = 1; 
%interactive_mode.rot_axis_tilt = 1;

tomo.rot_axis.offset = 5.4 * 2 / raw_bin; % 7.1.20
scan_path = [raw_path '008_WT_stage23_1_1']; 

image_correlation.area_height = [1 600];
image_correlation.area_width = [0 1];
write.parfolder = 'image_corr_roi_0_1_1_600';
ADD

%% Rot axis checked
write.parfolder = '';
interactive_mode.rot_axis_pos = 1;

% scan_path = [raw_path '001_WT_stage18_1_1_spiral']; ADD
% scan_path = [raw_path '001_WT_stage18_1_1_stepscan']; ADD
% scan_path = [raw_path '004_WT_stage23_1_1_spiral']; ADD

% vertical movement: ~5-8 pixels binned 2x
tomo.rot_axis.offset = 5.4 * 2 / raw_bin; % 7.1.20
scan_path = [raw_path '008_WT_stage23_1_1']; ADD

% vertical movement
tomo.rot_axis.offset = -3.95  * 2 / raw_bin; % 8.1.20
scan_path = [raw_path '009_WT_stage23_1_2']; ADD

% vertical movement
tomo.rot_axis.offset = 3.8  * 1 / raw_bin; % 8.1.20
scan_path = [raw_path '010_WT_stage22_1_1']; ADD

% vertical movement
tomo.rot_axis.offset = 4.8 * 2 / raw_bin; % 8.1.20
scan_path = [raw_path '011_WT_stage21_1_1']; ADD

% vertical movement
tomo.rot_axis.offset = 5.25 * 2 / raw_bin; % 8.1.20
scan_path = [raw_path '012_WT_stage22_2_1']; ADD

% vertical movment, 2-5 pixel unbinned
tomo.rot_axis.offset = 5.1 * 2 / raw_bin;
scan_path = [raw_path '013_Mohhmr_stage19_1_1']; ADD

% beam a bit more stable

% Poor quality, very weak contrast, flat field correction OK
% less vertical movment, but strong lateral movement or tilt or deformation
tomo.rot_axis.offset = 4.6 * 2 / raw_bin; % 8.1.20, more a guess, no peaks in metrics
scan_path = [raw_path '014_Mohmmr_stage23_1_1']; ADD

% scan_path = [raw_path '015_WT_stage18_1_1_spiral']; ADD
% scan_path = [raw_path '016_WT_stage18_1_1_spiral']; ADD
% scan_path = [raw_path '017_WT_stage18_1_1_spiral_360']; ADD
% scan_path = [raw_path '018_Mohmmr_stage18_1_1_spiral']; ADD
% scan_path = [raw_path '019_Mohmmr_stage18_2_1']; ADD

% vertical movment of bottom about 2-3 2x binned pixels
% breathing/deformation or lateral movement about 1-2 2x binnex pixels
%tomo.rot_axis.offset = 8.85 * 2 / raw_bin; % z = 0.2 % 10.1.20
tomo.rot_axis.offset = 8.25 * 2 / raw_bin; % z = 0.5
%tomo.rot_axis.offset = 7.5 * 2 / raw_bin; % z = 0.8
scan_path = [raw_path '020_Mohmmr_stage19_1_1']; ADD

% massive gas bubbles
% vertical movment about 1 2x binned pixels
image_correlation.area_height = [2000 3800];
image_correlation.area_width = [-1000 -201];
tomo.rot_axis.offset =  3.7 * 4 / raw_bin;
scan_path = [raw_path '021_Mohmmr_stage19_2_1']; ADD

% gas bubbles
% 3-4 2x pixels vertical movement
tomo.rot_axis.offset =  5.0 * 4 / raw_bin;
scan_path = [raw_path '022_WT_stage19_1_1']; ADD

% 1-2 2x pixels vertical movment
tomo.rot_axis.offset =  7.7 * 2 / raw_bin; % 15.01.20
phase_retrieval.apply = 1;
scan_path = [raw_path '023_WT_stage20_1_1']; ADD

% 7-9 2x pixels vertical movement
tomo.rot_axis.offset = 3.6  * 4 / raw_bin;
scan_path = [raw_path '024_WT_stage21_1_1']; ADD

% 3-5 2x pixels vertical movement
tomo.rot_axis.offset =  3.6 * 4 / raw_bin;
scan_path = [raw_path '025_WT_stage23_1_1']; ADD

% 1-3 2x pixels vertical movment, breathing
tomo.rot_axis.offset = 3.8  * 4 / raw_bin;
scan_path = [raw_path '026_Mohmmr_stage20_1_1']; ADD

% 1 2x pixels vertical movement, breathing
tomo.rot_axis.offset = 3.6  * 4 / raw_bin;
scan_path = [raw_path '027_Mohmmr_stage20_2_1']; ADD

% 5-7 2x pixels vertical movement, breathing
tomo.rot_axis.offset =  3.2 * 4 / raw_bin;
scan_path = [raw_path '028_Mohmmr_stage21_1_1']; ADD

% 1-4 2x pixels vertical movement
tomo.rot_axis.offset =  3.8 * 4 / raw_bin;
scan_path = [raw_path '029_Mohmmr_stage22_1_1']; ADD

% 5 2x pixel vertical movement
tomo.rot_axis.offset = 3.6 * 4 / raw_bin;
scan_path = [raw_path '030_Mohmmr_stage23_1_1']; ADD

% 5 2x pixel vertical movement
tomo.rot_axis.offset = 3.75 * 4 / raw_bin;
scan_path = [raw_path '031_Mohmmr_stage23_2_1']; ADD

% 2-3 2x pixel vertical movement, breathing, growing gas bubbles
tomo.rot_axis.offset = 4.0 * 4 / raw_bin;
scan_path = [raw_path '032_WT_stage18_1_1']; ADD

% 1-3 pixel 2x vertical movement, breathing
scan_path = [raw_path '033_Mohmmr_stage24_1_1']; ADD

% Only 1000 projections?
% 2-5 2x pixel vertical movement
scan_path = [raw_path '034_Mohmmr_stage19_3_1_IOD']; ADD

% 2-4 2x pixel vertical movement, breathing rather small compared to other
scan_path = [raw_path '035_Mohmmr_stage19_4_1_PTA']; ADD


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p05_reco_loop( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
