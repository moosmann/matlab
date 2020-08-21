function p05_reco_loop_hnee_tension_wood_000( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
% Created on 02-Jul-2018 by moosmanj

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
%%% SCAN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
read_flatcor = 0; % read preprocessed flatfield-corrected projections. CHECK if negative log has to be taken!
read_flatcor_path = ''; % subfolder of 'flat_corrected' containing projections
read_sino = 0; % read preprocessed sinograms. CHECK if negative log has to be taken!
read_sino_folder = ''; % subfolder to scan path
energy = []; % in eV! if empty: read from log file
sample_detector_distance = []; % in m. if empty: read from log file
eff_pixel_size = []; % in m. if empty: read from log lfile. effective pixel size =  detector pixel size / magnification
pixel_scaling = 1; % to account for beam divergence if pixel size was determined (via MTF) at the wrong distance
%%% PREPROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_roi = -1;[]; % vertical and/or horizontal ROI; (1,1) coordinate = top left pixel; supports absolute, relative, negative, and mixed indexing.
% []: use full image;
% [y0 y1]: vertical ROI, skips first raw_roi(1)-1 lines, reads until raw_roi(2); if raw_roi(2) < 0 reads until end - |raw_roi(2)|; relative indexing similar.
% [y0 y1 x0 x1]: vertical + horzontal ROI, each ROI as above
% if -1: auto roi, selects vertical ROI automatically. Use only for DCM. Not working for *.raw data where images are flipped and DMM data.
% if < -1: Threshold is set as min(proj(:,:,[1 end])) + abs(raw_roi)*median(dark(:)). raw_roi=-1 defaults to min(proj(:,:,[1 end])) + 4*median(dark(:))
raw_bin = 2; % projection binning factor: integer
im_trafo = '' ;%'rot90(im,-1)'; % string to be evaluated after reading data in the case the image is flipped/rotated/etc due to changes at the beamline, e.g. 'rot90(im)'
par.crop_at_rot_axis = 0; % for recos of scans with excentric rotation axis but WITHOUT projection stitching
par.stitch_projections = 0; % for 2 pi scans: stitch projection at rotation axis position. Recommended with phase retrieval to reduce artefacts. Standard absorption contrast data should work well without stitching. Subpixel stitching not supported (non-integer rotation axis position is rounded, less/no binning before reconstruction can be used to improve precision).
par.stitch_method = 'sine'; 'step';'linear'; %  ! CHECK correlation area !
% 'step' : no interpolation, use step function
% 'linear' : linear interpolation of overlap region
% 'sine' : sinusoidal interpolation of overlap region
proj_range = []; % range of projections to be used (from all found, if empty or 1: all, if scalar: stride, if range: start:incr:end
ref_range = []; % range of flat fields to be used (from all found), if empty or 1: all. if scalar: stride, if range: start:incr:end
pixel_filter_threshold_dark = [0.01 0.005]; % Dark fields: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_threshold_flat = [0.01 0.005]; % Flat fields: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_threshold_proj = [0.01 0.005]; % Raw projection: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_radius = [3 3]; % Increase only if blobs of zeros or other artefacts are expected. Can increase processing time heavily.
ring_current_normalization = 1; % normalize flat fields and projections by ring current
image_correlation.method = 'ssim';'entropy';'ssim-ml';'ssim-g';'std';'cov';'corr';'diff1-l1';'diff1-l2';'diff2-l1';'diff2-l2';'cross-entropy-12';'cross-entropy-21';'cross-entropy-x';
% Correlation of projections and flat fields. Essential for DCM data. Typically improves reconstruction quality of DMM data, too.
% Available methods ('ssim-ml'/'entropy' usually work best):
% 'none' : no correlation, uses median flat, for fast recos
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
phase_retrieval.reg_par = 0.5; % regularization parameter. larger values tend to blurrier images. smaller values tend to original data.
phase_retrieval.bin_filt = 0.1; % threshold for quasiparticle retrieval 'qp', 'qp2'
phase_retrieval.cutoff_frequ = 2 * pi; % in radian. frequency cutoff in Fourier space for 'qpcut' phase retrieval
phase_retrieval.padding = 1; % padding of intensities before phase retrieval, 0: no padding
%%% TOMOGRAPHY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tomo.run = 1; % run tomographic reconstruction
tomo.run_interactive_mode = 1; % if tomo.run = 0, use to determine rot axis positions without processing the full tomogram;
tomo.reco_mode = '3D'; 'slice'; % slice-wise or full 3D backprojection. 'slice': volume must be centered at origin & no support of rotation axis tilt, reco binning, save compressed
tomo.vol_size = []; %[-1 1 -1 1 -0.5 0.5];% 6-component vector [xmin xmax ymin ymax zmin zmax], for excentric rot axis pos / extended FoV;. if empty, volume is centerd within tomo.vol_shape. unit voxel size is assumed. if smaller than 10 values are interpreted as relative size w.r.t. the detector size. Take care bout minus signs! Note that if empty vol_size is dependent on the rotation axis position.
tomo.vol_shape = []; %[1 1 1] shape (# voxels) of reconstruction volume. used for excentric rot axis pos. if empty, inferred from 'tomo.vol_size'. in absolute numbers of voxels or in relative number w.r.t. the default volume which is given by the detector width and height.
tomo.rot_angle.full_range = []; % in radians. if []: full angle of rotation including additional increment, or array of angles. if empty full rotation angles is determined automatically to pi or 2 pi
tomo.rot_angle.offset = pi; % global rotation of reconstructed volume
tomo.rot_axis.offset = []; 
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
write.parfolder = 'for_russia';% parent folder to 'reco', 'sino', 'phase', and 'flat_corrected'
write.subfolder.flatcor = ''; % subfolder in 'flat_corrected'
write.subfolder.phase_map = ''; % subfolder in 'phase_map'
write.subfolder.sino = ''; % subfolder in 'sino'
write.subfolder.reco = ''; % subfolder in 'reco'
write.flatcor = 0; % save preprocessed flat corrected projections
write.phase_map = 0; % save phase maps (if phase retrieval is not 0)
write.sino = 1; % save sinograms (after preprocessing & before FBP filtering and phase retrieval)
write.phase_sino = 1; % save sinograms of phase maps
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
% dynamic range is compressed s.t. new dynamic range assumes
% 'outlier' : [LOW, HIGH] = write.compression.parameter, eg. [0.01 0.03], outlier given in percent, if scalear LOW = HIGH.
% 'full' : full dynamic range is used
% 'threshold' : [LOW HIGH] = write.compression.parameter, eg. [-0.01 1]
% 'std' : NUM = write.compression.parameter, mean +/- NUM*std, dynamic range is rescaled to within -/+ NUM standard deviations around the mean value
% 'histo' : [LOW HIGH] = write.compression.parameter (100*LOW)% and (100*HIGH)% of the original histogram, e.g. [0.02 0.02]
%%% INTERACTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.visual_output = 1; % show images and plots during reconstruction
interactive_mode.rot_axis_pos = 0; % reconstruct slices with dif+ferent rotation axis offsets
interactive_mode.rot_axis_pos_default_search_range = []; % if empty: asks for search range when entering interactive mode
interactive_mode.rot_axis_tilt = 0; % reconstruct slices with different offset AND tilts of the rotation axis
interactive_mode.rot_axis_tilt_default_search_range = []; % if empty: asks for search range when entering interactive mode
interactive_mode.lamino = 0; % find laminography tilt instead camera rotation
interactive_mode.fixed_other_tilt = 0; % fixed other tilt
interactive_mode.slice_number = 0.5; % default slice number. if in [0,1): relative, if in (1, N]: absolute
interactive_mode.phase_retrieval = 1; % Interactive retrieval to determine regularization parameter
interactive_mode.phase_retrieval_default_search_range = []; % if empty: asks for search range when entering interactive mode, otherwise directly start with given search range
%%% HARDWARE / SOFTWARE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.use_cluster = 1; % if available: on MAXWELL nodes disp/nova/wga/wgs cluster computation can be used. Recommended only for large data sets since parpool creation and data transfer implies a lot of overhead.
par.poolsize = 0.6; % number of workers used in a local parallel pool. if 0: use current config. if >= 1: absolute number. if 0 < poolsize < 1: relative amount of all cores to be used. if SLURM scheduling is available, a default number of workers is used.
par.poolsize_gpu_limit_factor = 0.7; % Relative amount of GPU memory used for preprocessing during parloop. High values speed up Proprocessing, but increases out-of-memory failure
tomo.astra_link_data = 1; % ASTRA data objects become references to Matlab arrays. Reduces memory issues.
tomo.astra_gpu_index = []; % GPU Device index to use, Matlab notation: index starts from 1. default: [], uses all
%%% EXPERIMENTAL OR NOT YET IMPLEMENTED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write.uint8_segmented = 0; % experimental: threshold segmentaion for histograms with 2 distinct peaks: __/\_/\__

% raw_roi = [301 3300];
% %%% PHASE RETRIEVAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% phase_retrieval.apply = 1; 
% phase_retrieval.apply_before = 1; 
% phase_retrieval.post_binning_factor = 1;
% phase_retrieval.method = 'tie';'qpcut'; 
% phase_retrieval.reg_par = 1.5; 
% phase_retrieval.bin_filt = 0.15;
% phase_retrieval.cutoff_frequ = 2 * pi;
% phase_retrieval.padding = 1;
% %%% TOMOGRAPHY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tomo.run = 1; 
% tomo.run_interactive_mode = 0;
% tomo.reco_mode = '3D'; 
% tomo.vol_size = [-0.5 0.5 -0.5 0.5 -0.5 0.5];
% tomo.vol_shape = [];
% tomo.rot_axis.offset = [];
% tomo.rot_axis.position = []; 
% tomo.rot_angle.full_range = []; 
% tomo.rot_angle.offset = pi / 2; 
% tomo.rot_axis.tilt = 0; 
% tomo.rot_axis.corr_area1 = [];
% tomo.rot_axis.corr_area2 = [];
% tomo.rot_axis.corr_gradient = 0;
% tomo.fbp_filter.type = 'linear';
% tomo.fbp_filter.freq_cutoff = 1;
% tomo.fbp_filter.padding = 1; 
% tomo.fbp_filter.padding_method = 'symmetric';
% tomo.butterworth_filter.apply = 1; % use butterworth filter in addition to FBP filter
% tomo.butterworth_filter.order = 1;
% tomo.butterworth_filter.frequ_cutoff = 0.95;
% tomo.astra_pixel_size = 1; 
% tomo.take_neg_log = []; 
% tomo.algorithm = 'fbp';'sirt'; 'cgls';
% tomo.iterations = 40; % for 'sirt' or 'cgls'.
% tomo.sirt.MinConstraint = []; 
% tomo.sirt.MaxConstraint = []; 

% Set default. Defines default parameter set to be used.
SET_DEFAULT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETER / DATA SETS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_path = '/asap3/petra3/gpfs/p05/2018/data/11004450/raw/';

%image_correlation.method = 'none';
%raw_roi = [2301 2600];
%write.reco = 0;
%tomo.run = 0;

%scan_path = [raw_path 'hnee01_kie_sh_ts_000N']; ADD
%scan_path = [raw_path 'hnee02_kie_sh_ts_000N']; ADD
%scan_path = [raw_path 'hnee03_kie_sh_ts_000N']; ADD
%scan_path = [raw_path 'hnee04_kie_sh_ts_000N']; ADD
%scan_path = [raw_path 'hnee05_kie_sh_ts_000N']; ADD
%scan_path = [raw_path 'hnee06_kie_sh_ts_000N']; ADD % some mismatch
%raw_roi = [301 3300];
tomo.rot_axis.offset = 2 * -784 / raw_bin;
scan_path = [raw_path 'hnee07_kie_sh_ts_000N']; ADD

%raw_roi = [301 3300];
tomo.rot_axis.offset = 2 * -783.5 / raw_bin;
scan_path = [raw_path 'hnee08_kie_sh_ts_000']; ADD
scan_path = [raw_path 'hnee09_kie_sh_ts_001']; ADD
scan_path = [raw_path 'hnee09_kie_sh_ts_002']; ADD
scan_path = [raw_path 'hnee09_kie_sh_ts_003']; ADD
scan_path = [raw_path 'hnee09_kie_sh_ts_004']; ADD
scan_path = [raw_path 'hnee09_kie_sh_ts_005']; ADD
scan_path = [raw_path 'hnee09_kie_sh_ts_006']; ADD
scan_path = [raw_path 'hnee09_kie_sh_ts_007']; ADD
scan_path = [raw_path 'hnee09_kie_sh_ts_008']; ADD
scan_path = [raw_path 'hnee09_kie_sh_ts_009']; ADD
scan_path = [raw_path 'hnee09_kie_sh_ts_010']; ADD

scan_path = [raw_path 'hnee10_kie_fh_ts_0p5N']; ADD
scan_path = [raw_path 'hnee11_kie_fh_ts_0p5N']; ADD
scan_path = [raw_path 'hnee12_kie_fh_ts_0p5N']; ADD
scan_path = [raw_path 'hnee13_kie_fh_ts_0p5N']; ADD
scan_path = [raw_path 'hnee14_kie_fh_ts_0p5N']; ADD

% not working
% scan_path = [raw_path 'hnee15_kie_fh_ts_000']; ADD
%scan_path = [raw_path 'hnee16_kie_fh_ts_000']; ADD
%scan_path = [raw_path 'hnee17_kie_fh_ts_000']; ADD

tomo.rot_axis.offset = 2 * -10.0 / raw_bin;
scan_path = [raw_path 'hnee18_pappel_tensionWood_000']; ADD
tomo.rot_axis.offset = 2 * -10.75 / raw_bin;
scan_path = [raw_path 'hnee18_pappel_tensionWood_001']; ADD
tomo.rot_axis.offset = 2 * -10.75 / raw_bin;
scan_path = [raw_path 'hnee18_pappel_tensionWood_002']; ADD
tomo.rot_axis.offset = 2 * -10.85 / raw_bin;
scan_path = [raw_path 'hnee18_pappel_tensionWood_003']; ADD
tomo.rot_axis.offset = 2 * -11.0 / raw_bin;
scan_path = [raw_path 'hnee18_pappel_tensionWood_004']; ADD
scan_path = [raw_path 'hnee18_pappel_tensionWood_005']; ADD
tomo.rot_axis.offset = 2 * -10.75 / raw_bin;
scan_path = [raw_path 'hnee18_pappel_tensionWood_006']; ADD
tomo.rot_axis.offset = 2 * -11.0 / raw_bin;
scan_path = [raw_path 'hnee18_pappel_tensionWood_007']; ADD
tomo.rot_axis.offset = 2 * -11.25 / raw_bin;
scan_path = [raw_path 'hnee18_pappel_tensionWood_008']; ADD

tomo.rot_axis.offset = 2 * -10.85 / raw_bin;
scan_path = [raw_path 'hnee19_pappel_oppositeWood_000']; ADD
tomo.rot_axis.offset = 2 * -10.3 / raw_bin;
scan_path = [raw_path 'hnee19_pappel_oppositeWood_001']; ADD
tomo.rot_axis.offset = 2 * -10.35 / raw_bin;
scan_path = [raw_path 'hnee19_pappel_oppositeWood_002']; ADD
tomo.rot_axis.offset = 2 * -10.35 / raw_bin;
scan_path = [raw_path 'hnee19_pappel_oppositeWood_003']; ADD 
% bad reco
scan_path = [raw_path 'hnee19_pappel_oppositeWood_004']; ADD
tomo.rot_axis.offset = 2 * -10.6 / raw_bin;
scan_path = [raw_path 'hnee19_pappel_oppositeWood_005']; ADD
tomo.rot_axis.offset = 2 * -10.7 / raw_bin;
scan_path = [raw_path 'hnee19_pappel_oppositeWood_006']; ADD
tomo.rot_axis.offset = 2 * -10.75 / raw_bin;
scan_path = [raw_path 'hnee19_pappel_oppositeWood_007']; ADD
tomo.rot_axis.offset = 2 * -10.75 / raw_bin;
scan_path = [raw_path 'hnee19_pappel_oppositeWood_008']; ADD

% Movement
tomo.rot_axis.offset = 2 * 8.5 / raw_bin;
scan_path = [raw_path 'hnee20_pappel_tensionWood_000']; ADD
tomo.rot_axis.offset = 2 * 8.5 / raw_bin;
scan_path = [raw_path 'hnee20_pappel_tensionWood_001']; ADD
tomo.rot_axis.offset = 2 * 8.25 / raw_bin;
scan_path = [raw_path 'hnee20_pappel_tensionWood_002']; ADD
tomo.rot_axis.offset = 2 * 8.0 / raw_bin;
scan_path = [raw_path 'hnee20_pappel_tensionWood_003']; ADD

tomo.rot_axis.offset = 2 * 7.5 / raw_bin;
scan_path = [raw_path 'hnee21_pappel_oppositeWood_000']; ADD
tomo.rot_axis.offset = 2 * 9.0 / raw_bin;
scan_path = [raw_path 'hnee21_pappel_oppositeWood_001']; ADD
tomo.rot_axis.offset = 2 * 8.85 / raw_bin;
scan_path = [raw_path 'hnee21_pappel_oppositeWood_002']; ADD
tomo.rot_axis.offset = 2 * 8.9 / raw_bin;
scan_path = [raw_path 'hnee21_pappel_oppositeWood_003']; ADD
tomo.rot_axis.offset = 2 * 8.9 / raw_bin;
scan_path = [raw_path 'hnee21_pappel_oppositeWood_004']; ADD
tomo.rot_axis.offset = 2 * 8.8 / raw_bin;
scan_path = [raw_path 'hnee21_pappel_oppositeWood_005']; ADD
tomo.rot_axis.offset = 2 * 9.05 / raw_bin;
scan_path = [raw_path 'hnee21_pappel_oppositeWood_006']; ADD

tomo.rot_axis.offset = 2 * 9.1 / raw_bin;
scan_path = [raw_path 'hnee22_pappel_oppositeWood_000']; ADD
tomo.rot_axis.offset = 2 * 9.1 / raw_bin;
scan_path = [raw_path 'hnee22_pappel_oppositeWood_001']; ADD
tomo.rot_axis.offset = 2 * 9.1 / raw_bin;
scan_path = [raw_path 'hnee22_pappel_oppositeWood_002']; ADD
tomo.rot_axis.offset = 2 * 10.3 / raw_bin;
scan_path = [raw_path 'hnee22_pappel_oppositeWood_003']; ADD
tomo.rot_axis.offset = 2 * 10.3 / raw_bin;
scan_path = [raw_path 'hnee22_pappel_oppositeWood_004']; ADD
tomo.rot_axis.offset = 2 * 10.39 / raw_bin;
scan_path = [raw_path 'hnee22_pappel_oppositeWood_005']; ADD

tomo.rot_axis.offset = 2 * 9.5 / raw_bin;
scan_path = [raw_path 'hnee23_pappel_tensionWood_000']; ADD
tomo.rot_axis.offset = 2 * 9.9 / raw_bin;
scan_path = [raw_path 'hnee23_pappel_tensionWood_001']; ADD
tomo.rot_axis.offset = 2 * 9.95 / raw_bin;
scan_path = [raw_path 'hnee23_pappel_tensionWood_002']; ADD
tomo.rot_axis.offset = 2 * 10.0 / raw_bin;
scan_path = [raw_path 'hnee23_pappel_tensionWood_003']; ADD
tomo.rot_axis.offset = 2 * 10.05 / raw_bin;
scan_path = [raw_path 'hnee23_pappel_tensionWood_004']; ADD
% No top up beam current drops from 92 to 82
tomo.rot_axis.offset = 2 * 10.10 / raw_bin;
scan_path = [raw_path 'hnee23_pappel_tensionWood_005']; ADD
% No top up beam current drops from 72? to 
tomo.rot_axis.offset = 2 * 10.05 / raw_bin;
scan_path = [raw_path 'hnee23_pappel_tensionWood_006']; ADD
% No top up beam current drops from 66 to 60
tomo.rot_axis.offset = 2 * 10.0 / raw_bin;
scan_path = [raw_path 'hnee23_pappel_tensionWood_007']; ADD
% No top up beam current drops from 58 to 55
tomo.rot_axis.offset = 2 * 10.0 / raw_bin;
scan_path = [raw_path 'hnee23_pappel_tensionWood_008']; ADD
% No top up beam current drops from 52 to 48
tomo.rot_axis.offset = 2 * 10.0 / raw_bin;
scan_path = [raw_path 'hnee23_pappel_tensionWood_009']; ADD
% No top up beam current drops from 47 to 44
tomo.rot_axis.offset = 2 * 10.0 / raw_bin;
scan_path = [raw_path 'hnee23_pappel_tensionWood_010']; ADD
% No top up beam current drops from 42 to 40
tomo.rot_axis.offset = 2 * 9.9 / raw_bin;
scan_path = [raw_path 'hnee23_pappel_tensionWood_011']; ADD
% No top up beam current drops from 38 to 36
tomo.rot_axis.offset = 2 * 9.95 / raw_bin;
scan_path = [raw_path 'hnee23_pappel_tensionWood_012']; ADD
% No top up beam current drops from 35 to 33
%% CHECK
tomo.rot_axis.offset = 2 * 10 / raw_bin;
scan_path = [raw_path 'hnee23_pappel_tensionWood_013']; ADD

%%%% TEST SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tomo.rot_axis.offset = 2 * 9.95 / raw_bin;
scan_path = [raw_path 'hnee23_pappel_tensionWood_002']; ADD

tomo.algorithm = 'sirt';
tomo.reco_mode = 'slice';
tomo.rot_axis.offset = 2 * 9.5 / raw_bin;
write.subfolder.reco = 'sirt'; 
scan_path = [raw_path 'hnee23_pappel_tensionWood_000']; ADD


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p05_reco_loop( SUBSETS, RUN_RECO, PRINT_PARAMETERS)