function p05_reco_loop_nano_11014414( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
% Created on 21-Jun-2022 by moosmanj

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
%%% SCAN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%l%%%%%%%%%%%
par.ref_path = {}; % cell of strings. Additonal data sets to be included for the correlation of projections and reference images
par.read_flatcor = 0; % read preprocessed flatfield-corrected projections. CHECK if negative log has to be taken!
par.read_flatcor_path = ''; % absolute path containing flat-field corrected projections
par.read_flatcor_trafo = @(im) im; %fliplr( im ); % anonymous function applied to the image which is read e.g. @(x) rot90(x)
par.read_sino = 0; % read preprocessed sinograms. CHECK if negative log has to be taken!
par.read_sino_folder = ''; % subfolder to scan path
par.read_sino_trafo = @(x) (x);%rot90(x); % anonymous function applied to the image which is read e.g. @(x) rot90(x)
par.energy = []; % eV! if empty: read from log file (log file values can be ambiguous or even missing sometimes)
par.sample_detector_distance = []; % in m. if empty: read from log file
par.eff_pixel_size = []; % in m. if empty: read from log lfile. effective pixel size =  detector pixel size / magnification
par.pixel_scaling = []; % to account for mismatch of eff_pixel_size with, ONLY APPLIED BEFORE TOMOGRAPHIC RECONSTRUCTION, HAS TO BE CHANGED!
par.read_image_log = 0; % bool, default: 0. Read metadata from image log instead hdf5, if image log exists
%%% PREPROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.raw_bin = 2; % projection binning factor: integer
par.raw_roi = []; % vertical and/or horizontal ROI; (1,1) coordinate = top left pixel; supports absolute, relative, negative, and mixed indexing.
% []: use full image;
% [y0 y1]: vertical ROI, skips first raw_roi(1)-1 lines, reads until raw_roi(2); if raw_roi(2) < 0 reads until end - |raw_roi(2)|; relative indexing similar.
% [y0 y1 x0 x1]: vertical + horzontal ROI, each ROI as above
% if scalar = 0: auto roi, selects vertical ROI automatically. Use only for DCM. Not working for *.raw data where images are flipped and DMM data.
% if scalar < 1: Threshold is set as min(proj(:,:,[1 end])) + abs(raw_roi)*median(dark(:)). raw_roi=-1 defaults to min(proj(:,:,[1 end])) + 4*median(dark(:))
par.im_trafo = '';% 'rot90(im,1)'; % string to be evaluated after reading data in the case the image is flipped/rotated/etc due to changes at the beamline, e.g. 'rot90(im)'
par.filter_ref = @(x) (x);
par.filter_proj = @(x) (x);
% STITCHING/CROPPING only for scans without lateral movment. Legacy support
par.crop_at_rot_axis = 0; % for recos of scans with excentric rotation axis but WITHOUT projection stitching
par.stitch_projections = 0; % for 2 pi cans: stitch projection at rotation axis position. Recommended with phase retrieval to reduce artefacts. Standard absorption contrast data should work well without stitching. Subpixel stitching not supported (non-integer rotation axis position is rounded, less/no binning before reconstruction can be used to improve precision).
par.stitch_method = 'sine'; 'step';'linear'; %  ! CHECK correlation area !
% 'step' : no interpolation, use step function
% 'linear' : linear interpolation of overlap region
% 'sine' : sinusoidal interpolation of overlap region
par.proj_range = []; % range of projections to be used (from all found, if empty or 1: all, if scalar: stride, if range: start:incr:end
par.ref_range = [];%1:100;[]; % range of flat fields to be used (from all found), if empty or 1: all. if scalar: stride, if range: start:incr:end
par.wo_crop = 0; % Do not crop images to account for random lateral shift
par.virt_s_pos = 0; % Correct sample position in reconsructed volume if virtual sample position motors are used
pixel_filter_threshold_dark = [0.01 0.005]; % Dark fields: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_threshold_flat = [0.02 0.005]; % Flat fields: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_threshold_proj = [0.02 0.005]; % Raw projection: threshold parameter for hot/dark pixel filter, for details see 'FilterPixel'
pixel_filter_radius = [5 5]; % Increase only if blobs of zeros or other artefacts are expected. Can increase processing time heavily.
par.ring_current_normalization = 1; % normalize flat fields and projections by ring current
image_correlation.method = 'ssim-ml';'median';'entropy';'none';'ssim';'ssim-g';'std';'cov';'corr';'diff1-l1';'diff1-l2';'diff2-l1';'diff2-l2';'cross-entropy-12';'cross-entropy-21';'cross-entropy-x';
image_correlation.force_calc = 0; % bool. force compuation of correlation even though a (previously computed) corrlation matrix exists
image_correlation.num_flats = 11; % number of best maching flat fields used for correction
image_correlation.area_width = [1 100];%[-100 1];% correlation area: index vector or relative/absolute position of [first pix, last pix], negative indexing is supported
image_correlation.area_height = [0.2 0.8]; % correlation area: index vector or relative/absolute position of [first pix, last pix]
image_correlation.filter = 1; % bool, filter ROI before correlation
image_correlation.filter_type = 'median'; % string. correlation ROI filter type, currently only 'median' is implemnted
image_correlation.filter_parameter = {[3 3], 'symmetric'}; % cell. filter paramaters to be parsed with {:}
ring_filter.apply_before_stitching = 1; % ! Consider when phase retrieval is applied !
ring_filter.method = 'jm'; 'wavelet-fft';
ring_filter.waveletfft_dec_levels = 1:6; % decomposition levels for 'wavelet-fft'
ring_filter.waveletfft_wname = 'db7';'db25';'db30'; % wavelet type, see 'FilterStripesCombinedWaveletFFT' or 'waveinfo'
ring_filter.waveletfft_sigma = 3; %  suppression factor for 'wavelet-fft'
ring_filter.jm_median_width = 11; % multiple widths are applied consecutively, eg [3 11 21 31 39];
par.strong_abs_thresh = 1; % if 1: does nothing, if < 1: flat-corrected values below threshold are set to one. Try with algebratic reco techniques.
par.norm_sino = 0; % not recommended, can introduce severe artifacts, but sometimes improves quality
%%% PHASE RETRIEVAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phase_retrieval.apply_before = 0; % before stitching, interactive mode, etc. For phase-contrast data with an excentric rotation axis phase retrieval should be done afterwards. To find the rotataion axis position use this option in a first run, and then turn it of afterwards.
phase_retrieval.post_binning_factor = 1; % Binning factor after phase retrieval, but before tomographic reconstruction
phase_retrieval.method = 'tie';'tieNLO_Schwinger';'dpc';'tie';'qp';'qpcut'; %'qp' 'ctf' 'tie' 'qp2' 'qpcut'
% Interactive phase retrieval not supported for method 'tieNLO_Schwinger'
phase_retrieval.reg_par = 2; % regularization parameter. larger values tend to blurrier images. smaller values tend to original data.
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
tomo.vol_size = [];[-1 1 -1 1 -0.5 0.5];% 6-component vector [xmin xmax ymin ymax zmin zmax]for excentric rot axis pos / extended FoV;. if empty, volume is centerd within tomo.vol_shape. unit voxel size is assumed. if smaller than 10 values are interpreted as relative size w.r.t. the detector size. Take care bout minus signs! Note that if empty vol_size is dependent on the rotation axis position.
% [ slice left,  slice right,  slice top,  slice bottom, bottom slice, top slice], 2D Matlab notation: (0,0)=topleft:
tomo.vol_shape = []; %[1 1 1] shape (# voxels) of reconstruction volume. used for excentric rot axis pos. if empty, inferred from 'tomo.vol_size'. in absolute numbers of voxels or in relative number w.r.t. the default volume which is given by the detector width and height.
tomo.rot_angle_full_range = [];% 2 * pi ;[]; % in radians. if []: full angle of rotation including additional increment, or array of angles. if empty full rotation angles is determined automatically to pi or 2 pi
tomo.rot_angle_offset = 0; % global rotation of reconstructed volume
tomo.interpolate_missing_angles = 0; % limited or missing angle tomography
tomo.rot_axis_offset = [] / 2 * par.raw_bin; % rotation axis offset w.r.t to the image center. Assuming the rotation axis position to be centered in the FOV for standard scan, the offset should be close to zero.
tomo.rot_axis_position = []; % if empty use automatic computation. EITHER OFFSET OR posITION MUST BE EMPTY. YOU MUST NOT USE BOTH!
tomo.rot_axis_offset_shift = []; %[]; % absolute lateral movement in pixels during fly-shift-scan, overwrite lateral shift read out from hdf5 log
tomo.vert_shift = [];
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
tomo.rot_axis_search_auto = 0; % find extrema of metric within search range
tomo.rot_axis_search_range = []; % search reach for automatic determination of the rotation axis offset, overwrite interactive result if not empty
tomo.rot_axis_search_metric = 'neg'; % string: 'neg','entropy','iso-grad','laplacian','entropy-ML','abs'. Metric to find rotation axis offset
tomo.rot_axis_search_extrema = 'max'; % string: 'min'/'max'. chose min or maximum position
tomo.rot_axis_search_fit = 1; % bool: fit calculated metrics and find extrema, otherwise use extrema from search range
%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
write.path = ''; %'/gpfs/petra3/scratch/moosmanj'; % absolute path were output data will be stored. !!overwrites the write.to_scratch flag. if empty uses the beamtime directory and either 'processed' or 'scratch_cc'
write.deleteFiles = 0; % delete files already existing in output folders. Useful if number or names of files differ when reprocessing.
write.beamtimeID = ''; % string (regexp),typically beamtime ID, mandatory if 'write.deleteFiles' is true (safety check)
write.scan_name_appendix = ''; % appendix to the output folder name which defaults to the scan name
write.parfolder = ''; % parent folder to 'reco', 'sino', 'phase', and 'flat_corrected'
write.subfolder_flatcor = ''; % subfolder in 'flat_corrected'
write.subfolder_phase_map = ''; % subfolder in 'phase_map'
write.subfolder_sino = ''; % subfolder in 'sino'
write.subfolder_reco = ''; % subfolder in 'reco'
write.flatcor = 1; % save preprocessed flat corrected projections
write.phase_map = 0; % save phase maps (if phase retrieval is not 0)
write.sino = 0; % save sinograms (after preprocessing & before FBP filtering and phase retrieval)
write.phase_sino = 0; % save sinograms of phase maps
write.reco = 1; % save reconstructed slices (if tomo.run=1)
write.float = 1; % single precision (32-bit float) tiff
write.uint16 = 0; % save 16bit unsigned integer tiff using 'write.compression_method'
write.uint8 = 0; % save binned 8bit unsigned integer tiff using 'write.compression_method'
% Optionally save binned reconstructions, only works in '3D' reco_mode
write.float_binned = 1; % save binned single precision (32-bit float) tiff
write.uint16_binned = 0; % save binned 16bit unsigned integer tiff using 'write.compression_method'
write.uint8_binned = 0; % save binned 8bit unsigned integer tiff using 'wwrite.compression_method'
write.reco_binning_factor = 2; % IF BINNED VOLUMES ARE SAVED: binning factor of reconstructed volume
write.compression_method = 'outlier';'histo';'full'; 'std'; 'threshold'; % method to compression dynamic range into [0, 1]
write.compression_parameter = [0.02 0.02]; % compression-method specific parameters
write.uint8_segmented = 0; % experimental: threshold segmentaion for histograms with 2 distinct peaks: __/\_/\__
%%% INTERACTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.visual_output = 1; % show images and plots during reconstruction
par.skip_gpu_info = 1;
interactive_mode.rot_axis_pos = 1; % reconstruct slices with dif+ferent rotation axis offsets
interactive_mode.rot_axis_pos_default_search_range = []; % if empty: asks for search range when entering interactive mode
interactive_mode.rot_axis_tilt = 0; % reconstruct slices with different offset AND tilts of the rotation axis
interactive_mode.rot_axis_tilt_default_search_range = []; % if empty: asks for search range when entering interactive mode
interactive_mode.lamino = 1; % find laminography tilt instead camera tilt
interactive_mode.angles = 0; % reconstruct slices with different scalings of angles
interactive_mode.angle_scaling_default_search_range = []; % if empty: use a variaton of -/+5 * (angle increment / maximum angle)
interactive_mode.slice_number = 0.5; % default slice number. if in [0,1): relative, if in (1, N]: absolute
interactive_mode.phase_retrieval_default_search_range = []; % if empty: asks for search range when entering interactive mode, otherwise directly start with given search range
%%% HARDWARE / SOFTWARE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.use_cluster = 0; % if available: on MAXWELL nodes disp/nova/wga/wgs cluster computation can be used. Recommended only for large data sets since parpool creation and data transfer implies a lot of overhead.
par.poolsize = 0.8; % number of workers used in a local parallel pool. if 0: use current config. if >= 1: absolute number. if 0 < poolsize < 1: relative amount of all cores to be used. if SLURM scheduling is available, a default number of workers is used.
par.poolsize_gpu_limit_factor = 0.8; % Relative amount of GPU memory used for preprocessing during parloop. High values speed up Proprocessing, but increases out-of-memory failure
tomo.astra_link_data = 1; % ASTRA data objects become references to Matlab arrays. Reduces memory issues.
tomo.astra_gpu_index = []; % GPU Device index to use, Matlab notation: index starts from 1. default: [], uses all
par.gpu_index = tomo.astra_gpu_index;
par.use_gpu_in_parfor = 1;
phase_retrieval.use_parpool = 1; % bool. Disable parpool when out-of-memory error occurs during phase retrieval.
par.window_state = 'minimized';'normal';'maximized';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set default. Defines default parameter set to be used.
SET_DEFAULT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETER / DATA SETS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_path = '/asap3/petra3/gpfs/p05/2022/data/11014414/raw/';
write.to_scratch = 0; 
ring_filter.apply = 0; 


par.scan_path = [raw_path 'imke_001_pocke']; ADD
par.scan_path = [raw_path 'imke_002_pocke']; ADD
par.scan_path = [raw_path 'imke_003_pocke']; ADD
write.to_scratch = 0;

par.read_filenames_from_disk = 1;
par.ring_current_normalization = 0;
phase_retrieval.apply = 0;
par.raw_roi = [1000, 4500];

interactive_mode.rot_axis_pos = 1;
interactive_mode.angles = 0; 

par.raw_bin = 2;
par.proj_range = 1:5501-82;
%eff_pix_size = 0.46;
%s = - s_stage_x / eff_pix_size * 1000;
s = s_stage_x;
tomo.rot_axis_offset_shift = s(par.proj_range);
a = 4 * pi * (0:5500)/5501;
par.crop_proj = 1;
tomo.rot_angle_full_range = a(par.proj_range);

par.scan_path = [raw_path 'imke_004_kaktus_golden1']; ADD

phase_retrieval.apply = 1;
interactive_mode.phase_retrieval = 1; 
phase_retrieval.method ='tie';'tieNLO_Schwinger';'dpc';'tie'; 'qp';'qpcut';

ADD
par.scan_path = [raw_path 'imke_005_kaktus_golden1']; ADD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p05_reco_loop( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
end

function s = s_stage_x()
    s = [  0.22528772      0.20247802      0.17145682      0.18787981      0.14271659
      0.21935720      0.19016078      0.16735107      0.14590995      0.19107316
      0.20384660      0.19335413      0.19107316      0.22848108      0.20932093
      0.17236921      0.15184047      0.15594622      0.17008824      0.17054443
      0.19837227      0.16598249      0.19654749      0.15640241      0.15594622
      0.23121824      0.20156563      0.18924839      0.14636614      0.20795235
      0.19198555      0.19472272      0.18970458      0.18742361      0.15685861
      0.19974085      0.22665630      0.16917585      0.21890100      0.17236921
      0.16050816      0.15731480      0.15092808      0.21570765      0.22620011
      0.15001570      0.16780727      0.19563511      0.18696742      0.21570765
      0.18423025      0.22984966      0.15092808      0.14408517      0.14773473
      0.22665630      0.14682234      0.17236921      0.14910331      0.21616384
      0.16963204      0.23213063      0.16826346      0.15594622      0.20567138
      0.16096435      0.20156563      0.15685861      0.17647495      0.14773473
      0.22620011      0.17966831      0.16643868      0.14317278      0.18240548
      0.15229667      0.20247802      0.16963204      0.18194928      0.22209436
      0.18058070      0.16005197      0.20612757      0.16050816      0.18651122
      0.19882846      0.14955950      0.19061697      0.19700369      0.16963204
      0.20977712      0.18696742      0.23121824      0.19837227      0.22665630
      0.22255056      0.16324532      0.22711250      0.22620011      0.17921212
      0.15412144      0.18423025      0.21433906      0.18742361      0.19609130
      0.18468645      0.15913958      0.14682234      0.16598249      0.14180420
      0.15320905      0.22391914      0.20475899      0.17875592      0.20293421
      0.21753242      0.17510637      0.19791608      0.19289794      0.20156563
      0.22026959      0.14819092      0.17373779      0.19517891      0.16643868
      0.22483153      0.18058070      0.15092808      0.20247802      0.22300675
      0.20977712      0.18377406      0.16005197      0.15047189      0.17556257
      0.23167444      0.16552630      0.14180420      0.20749615      0.18103689
      0.15640241      0.16370152      0.17100062      0.17875592      0.19016078
      0.16780727      0.18742361      0.22756869      0.22483153      0.20658376
      0.20703996      0.14864711      0.20384660      0.16507010      0.21981339
      0.22118197      0.17373779      0.20612757      0.22026959      0.15138428
      0.18833600      0.17510637      0.18514264      0.20019705      0.16050816
      0.15275286      0.15092808      0.19016078      0.17419398      0.19654749
      0.16324532      0.18742361      0.22209436      0.15640241      0.22984966
      0.19061697      0.19837227      0.20840854      0.20019705      0.23030586
      0.17145682      0.21114570      0.18377406      0.18194928      0.20567138
      0.19745988      0.20977712      0.18559884      0.19061697      0.17829973
      0.14454137      0.19472272      0.19654749      0.16324532      0.17738734
      0.15092808      0.17647495      0.22984966      0.20658376      0.19609130
      0.21433906      0.16689488      0.16917585      0.21342668      0.18651122
      0.21753242      0.19335413      0.15913958      0.16689488      0.15868338
      0.20977712      0.19974085      0.16096435      0.15959577      0.18468645
      0.21023332      0.21023332      0.22391914      0.16233294      0.17556257
      0.22300675      0.21662003      0.20156563      0.20886473      0.18514264
      0.19517891      0.21570765      0.22984966      0.18377406      0.18514264
      0.17236921      0.21890100      0.17465018      0.19700369      0.16050816
      0.14727853      0.22391914      0.17328159      0.22391914      0.22026959
      0.22711250      0.19700369      0.19517891      0.20658376      0.17784354
      0.20156563      0.15913958      0.21388287      0.22072578      0.19107316
      0.17829973      0.16415771      0.21935720      0.15640241      0.18058070
      0.17191301      0.22711250      0.17054443      0.18103689      0.23121824
      0.17465018      0.17601876      0.17601876      0.20521518      0.18103689
      0.14910331      0.16917585      0.20886473      0.18331786      0.19791608
      0.16963204      0.20430279      0.20293421      0.22528772      0.23076205
      0.22620011      0.22665630      0.17054443      0.23121824      0.16598249
      0.18970458      0.22893727      0.18696742      0.16643868      0.20065324
      0.14545376      0.15777100      0.17100062      0.20156563      0.19381033
      0.18240548      0.17556257      0.15229667      0.23030586      0.21388287
      0.14454137      0.21525145      0.19198555      0.14317278      0.18879219
      0.15047189      0.18559884      0.15138428      0.15868338      0.22939347
      0.14590995      0.18423025      0.19563511      0.20795235      0.22346295
      0.22209436      0.14682234      0.14271659      0.21935720      0.16598249
      0.15092808      0.20430279      0.19837227      0.15503383      0.18058070
      0.16096435      0.17191301      0.20293421      0.20293421      0.19244175
      0.22391914      0.15184047      0.18696742      0.20475899      0.18103689
      0.14727853      0.16780727      0.20977712      0.19609130      0.21753242
      0.22163817      0.23213063      0.22026959      0.19381033      0.18514264
      0.15549003      0.21753242      0.19837227      0.15685861      0.15777100
      0.17100062      0.18605503      0.20156563      0.14727853      0.15092808
      0.21798862      0.20156563      0.17738734      0.20065324      0.21570765
      0.22300675      0.18696742      0.20567138      0.15549003      0.20110943
      0.14180420      0.15275286      0.20840854      0.17738734      0.14499756
      0.21707623      0.18559884      0.17145682      0.22574392      0.19426652
      0.22026959      0.15366525      0.20475899      0.23167444      0.20475899
      0.16096435      0.15594622      0.20658376      0.18331786      0.22802489
      0.17647495      0.20932093      0.22574392      0.22163817      0.15457764
      0.15959577      0.14727853      0.18742361      0.19837227      0.20703996
      0.21068951      0.20795235      0.18240548      0.17875592      0.15594622
      0.21570765      0.17601876      0.14727853      0.17738734      0.20156563
      0.14408517      0.15731480      0.19791608      0.22711250      0.15412144
      0.20339040      0.14408517      0.20430279      0.15913958      0.23030586
      0.19061697      0.15275286      0.16826346      0.19289794      0.15913958
      0.18423025      0.18103689      0.17966831      0.16871965      0.16050816
      0.17738734      0.20293421      0.17465018      0.14819092      0.17008824
      0.20932093      0.15868338      0.15594622      0.17829973      0.15913958
      0.18879219      0.17784354      0.18787981      0.18423025      0.15777100
      0.16096435      0.16233294      0.17601876      0.21844481      0.18787981
      0.22163817      0.16735107      0.17328159      0.18240548      0.22620011
      0.18377406      0.20703996      0.20567138      0.19289794      0.14682234
      0.19928466      0.21798862      0.16187674      0.21981339      0.19974085
      0.16963204      0.20339040      0.16142055      0.18970458      0.20795235
      0.19882846      0.16370152      0.22528772      0.14590995      0.18331786
      0.19928466      0.21297048      0.18286167      0.17145682      0.18240548
      0.20795235      0.19654749      0.22528772      0.15184047      0.14271659
      0.17054443      0.20977712      0.17556257      0.20339040      0.20703996
      0.15549003      0.15366525      0.22209436      0.15868338      0.18058070
      0.19928466      0.18468645      0.18149309      0.15229667      0.16871965
      0.18833600      0.18696742      0.16415771      0.16233294      0.16826346
      0.20065324      0.14773473      0.18058070      0.18742361      0.18377406
      0.21342668      0.18924839      0.18331786      0.17784354      0.19061697
      0.16233294      0.14773473      0.20156563      0.15138428      0.20612757
      0.19335413      0.17966831      0.19517891      0.20521518      0.20886473
      0.15275286      0.18423025      0.19426652      0.16324532      0.20110943
      0.20703996      0.19837227      0.15412144      0.20293421      0.18605503
      0.19928466      0.15777100      0.16370152      0.20521518      0.16370152
      0.22255056      0.22437533      0.17328159      0.22574392      0.17054443
      0.20156563      0.18651122      0.18651122      0.21616384      0.15229667
      0.18058070      0.16598249      0.18012451      0.21342668      0.22620011
      0.18377406      0.19335413      0.21844481      0.14682234      0.17373779
      0.18559884      0.17373779      0.22620011      0.22437533      0.20795235
      0.19791608      0.15138428      0.20567138      0.15320905      0.21251429
      0.19426652      0.16963204      0.16050816      0.16005197      0.18331786
      0.18103689      0.14682234      0.22483153      0.15913958      0.17510637
      0.18149309      0.20293421      0.20567138      0.19472272      0.17601876
      0.14682234      0.17921212      0.22893727      0.16187674      0.15777100
      0.15959577      0.22437533      0.16552630      0.15229667      0.14180420
      0.15275286      0.18423025      0.19837227      0.15594622      0.22756869
      0.16826346      0.20293421      0.16507010      0.18331786      0.19381033
      0.15092808      0.18468645      0.14864711      0.19974085      0.21023332
      0.18879219      0.19745988      0.18286167      0.19381033      0.16050816
      0.14362898      0.18879219      0.22163817      0.18377406      0.15868338
      0.21662003      0.22528772      0.18833600      0.21205809      0.14819092
      0.18012451      0.17419398      0.17829973      0.20430279      0.20384660
      0.15959577      0.22118197      0.20567138      0.18879219      0.21981339
      0.21205809      0.15777100      0.15959577      0.14590995      0.16187674
      0.14955950      0.19061697      0.19198555      0.14819092      0.17966831
      0.14727853      0.21844481      0.18423025      0.21205809      0.20247802
      0.22848108      0.21388287      0.16415771      0.18012451      0.15184047
      0.21433906      0.21479526      0.16643868      0.17236921      0.20065324
      0.14180420      0.20156563      0.14454137      0.16278913      0.19609130
      0.14819092      0.14180420      0.21981339      0.15092808      0.23030586
      0.15822719      0.16507010      0.19882846      0.14819092      0.20110943
      0.16917585      0.14590995      0.19700369      0.14226040      0.15731480
      0.14271659      0.14636614      0.14454137      0.17191301      0.16780727
      0.20065324      0.17100062      0.19654749      0.21935720      0.17008824
      0.19472272      0.18286167      0.22391914      0.16233294      0.20293421
      0.18103689      0.16415771      0.21525145      0.16552630      0.17738734
      0.19198555      0.19745988      0.14317278      0.21616384      0.22574392
      0.16735107      0.15685861      0.15640241      0.19974085      0.18194928
      0.18514264      0.23167444      0.15685861      0.20567138      0.19472272
      0.16370152      0.20202182      0.19745988      0.22528772      0.20430279
      0.20795235      0.18787981      0.17236921      0.22300675      0.20612757
      0.17100062      0.22574392      0.18058070      0.18286167      0.16096435
      0.17465018      0.15913958      0.20247802      0.17008824      0.14682234
      0.17738734      0.17282540      0.21616384      0.20339040      0.18377406
      0.15913958      0.20430279      0.17191301      0.19700369      0.21205809
      0.14180420      0.18194928      0.16917585      0.22072578      0.23167444
      0.22483153      0.17829973      0.19700369      0.22118197      0.17556257
      0.14773473      0.19882846      0.23121824      0.14226040      0.18787981
      0.21342668      0.22483153      0.17875592      0.15457764      0.16142055
      0.16233294      0.16050816      0.17100062      0.20749615      0.17328159
      0.21798862      0.16963204      0.15503383      0.19198555      0.14454137
      0.17282540      0.17829973      0.21479526      0.20567138      0.15777100
      0.21251429      0.18149309      0.14819092      0.14955950      0.14545376
      0.15366525      0.20567138      0.16461391      0.18331786      0.19107316
      0.15001570      0.16005197      0.19882846      0.22300675      0.19517891
      0.19700369      0.20339040      0.18331786      0.17510637      0.19882846
      0.16461391      0.15822719      0.15229667      0.19426652      0.20612757
      0.14271659      0.19837227      0.17921212      0.15184047      0.21525145
      0.22346295      0.17829973      0.20475899      0.21342668      0.22528772
      0.14819092      0.22574392      0.14955950      0.19335413      0.17647495
      0.22391914      0.22209436      0.18286167      0.18194928      0.23030586
      0.21935720      0.16552630      0.16461391      0.15777100      0.19928466
      0.17100062      0.18742361      0.20703996      0.21023332      0.16643868
      0.19335413      0.15868338      0.19517891      0.15913958      0.18423025
      0.16871965      0.20019705      0.21616384      0.18696742      0.21890100
      0.22026959      0.15822719      0.19517891      0.15777100      0.20156563
      0.20886473      0.21707623      0.14773473      0.17328159      0.21433906
      0.22118197      0.15731480      0.20703996      0.14180420      0.15731480
      0.14955950      0.18286167      0.23030586      0.20612757      0.20703996
      0.19700369      0.21114570      0.20384660      0.22163817      0.18742361
      0.17921212      0.14864711      0.16096435      0.17647495      0.19061697
      0.20703996      0.17373779      0.20612757      0.14362898      0.15640241
      0.19381033      0.20567138      0.19335413      0.22163817      0.15503383
      0.18970458      0.21114570      0.23167444      0.17966831      0.17419398
      0.22802489      0.18696742      0.15959577      0.18286167      0.17191301
      0.17465018      0.15457764      0.21844481      0.16826346      0.14545376
      0.18377406      0.14636614      0.16370152      0.20521518      0.14590995
      0.21707623      0.21844481      0.16643868      0.15366525      0.19928466
      0.18194928      0.19654749      0.15092808      0.15457764      0.18194928
      0.18696742      0.18012451      0.15138428      0.23030586      0.17100062
      0.16735107      0.21023332      0.15731480      0.14864711      0.17100062
      0.23030586      0.20840854      0.20110943      0.22893727      0.19517891
      0.20612757      0.20840854      0.15959577      0.23121824      0.16370152
      0.15822719      0.15001570      0.16415771      0.20703996      0.22163817
      0.22711250      0.15913958      0.21205809      0.14408517      0.18377406
      0.16871965      0.21662003      0.21205809      0.21890100      0.17601876
      0.19563511      0.23121824      0.17465018      0.19198555      0.15001570
      0.22391914      0.23213063      0.17008824      0.17693115      0.19335413
      0.18103689      0.15685861      0.14636614      0.18377406      0.21068951
      0.17465018      0.15777100      0.14773473      0.18149309      0.17191301
      0.23213063      0.19882846      0.20932093      0.22118197      0.22346295
      0.22391914      0.15366525      0.14362898      0.22574392      0.22893727
      0.19198555      0.19016078      0.14499756      0.16871965      0.18468645
      0.14362898      0.21935720      0.15731480      0.19837227      0.19745988
      0.19517891      0.18605503      0.18286167      0.17556257      0.15275286
      0.21068951      0.23213063      0.20339040      0.15229667      0.20658376
      0.22802489      0.15549003      0.23076205      0.20339040      0.20384660
      0.16689488      0.15685861      0.18468645      0.23076205      0.16142055
      0.16552630      0.16643868      0.14362898      0.22528772      0.20840854
      0.21798862      0.22528772      0.18286167      0.16005197      0.14408517
      0.18514264      0.17008824      0.18833600      0.18286167      0.22209436
      0.20840854      0.20658376      0.21616384      0.19609130      0.14819092
      0.19837227      0.21798862      0.21570765      0.20293421      0.16826346
      0.20749615      0.14773473      0.16005197      0.19107316      0.14727853
      0.20886473      0.15275286      0.18468645      0.17054443      0.22528772
      0.23076205      0.21479526      0.18879219      0.21388287      0.18742361
      0.23030586      0.22574392      0.19152936      0.19882846      0.15184047
      0.19289794      0.17601876      0.17191301      0.22118197      0.15229667
      0.17282540      0.19654749      0.21160190      0.14773473      0.20612757
      0.16415771      0.18377406      0.20475899      0.20339040      0.17693115
      0.18651122      0.18377406      0.20156563      0.18058070      0.18194928
      0.19928466      0.19244175      0.17145682      0.21525145      0.15275286
      0.15001570      0.19426652      0.21068951      0.14819092      0.16917585
      0.16005197      0.19244175      0.14955950      0.20065324      0.17875592
      0.16370152      0.16507010      0.21570765      0.17738734      0.18742361
      0.21935720      0.22346295      0.21251429      0.15229667      0.21707623
      0.15685861      0.16461391      0.15959577      0.21662003      0.16233294
      0.16187674      0.15457764      0.19745988      0.23213063      0.20932093
      0.21068951      0.21205809      0.20202182      0.19700369      0.20339040
      0.21707623      0.18331786      0.20156563      0.17282540      0.14545376
      0.21297048      0.22437533      0.21981339      0.20110943      0.15184047
      0.17784354      0.20247802      0.20293421      0.17282540      0.20977712
      0.17601876      0.21753242      0.21890100      0.18605503      0.19016078
      0.17510637      0.21023332      0.17191301      0.18651122      0.20247802
      0.18423025      0.15640241      0.17465018      0.14819092      0.17054443
      0.18058070      0.18742361      0.15959577      0.14636614      0.16598249
      0.19974085      0.19700369      0.22255056      0.22163817      0.22528772
      0.15047189      0.17921212      0.19928466      0.14727853      0.20247802
      0.16917585      0.14682234      0.14955950      0.22893727      0.20247802
      0.19791608      0.16871965      0.22620011      0.15503383      0.20749615
      0.15138428      0.21753242      0.23121824      0.18286167      0.15184047
      0.20019705      0.20749615      0.20749615      0.20475899      0.19335413
      0.16643868      0.20886473      0.18970458      0.16050816      0.16598249
      0.22528772      0.17191301      0.19335413      0.14955950      0.22802489
      0.20339040      0.17465018      0.14454137      0.21068951      0.19426652
      0.18103689      0.21297048      0.22118197      0.16689488      0.20749615
      0.20110943      0.18696742      0.15184047      0.18605503      0.15092808
      0.21707623      0.19289794      0.19745988      0.19609130      0.23121824
      0.20110943      0.16598249      0.14910331      0.18103689      0.17236921
      0.18696742      0.16415771      0.22300675      0.15822719      0.18240548
      0.22026959      0.22163817      0.21342668      0.19837227      0.21433906
      0.20293421      0.15229667      0.20703996      0.17373779      0.18605503
      0.15503383      0.18696742      0.23030586      0.21890100      0.16096435
      0.22893727      0.20475899      0.17693115      0.20202182      0.20612757
      0.20339040      0.18468645      0.15320905      0.15868338      0.15366525
      0.23030586      0.16552630      0.20110943      0.16050816      0.14454137
      0.19061697      0.18149309      0.22163817      0.14180420      0.23167444
      0.16735107      0.17054443      0.21068951      0.19152936      0.21662003
      0.15366525      0.15959577      0.20521518      0.14819092      0.18924839
      0.18103689      0.18012451      0.21479526      0.14727853      0.14271659
      0.16552630      0.19472272      0.17784354      0.21160190      0.15457764
      0.16461391      0.22300675      0.14454137      0.17100062      0.16917585
      0.18468645      0.18559884      0.16643868      0.16735107      0.15913958
      0.14955950      0.20475899      0.18468645      0.14180420      0.23167444
      0.15457764      0.22848108      0.17236921      0.22756869      0.17693115
      0.18423025      0.20703996      0.18879219      0.19928466      0.15913958
      0.22939347      0.20977712      0.15229667      0.19016078      0.17738734
      0.15685861      0.21798862      0.16142055      0.16050816      0.15184047
      0.17966831      0.18423025      0.14590995      0.17145682      0.21707623
      0.15549003      0.14317278      0.18377406      0.14180420      0.23121824
      0.23121824      0.20521518      0.15457764      0.18651122      0.15184047
      0.15320905      0.19700369      0.14317278      0.21297048      0.20749615
      0.21205809      0.23167444      0.21205809      0.16963204      0.20110943
      0.18514264      0.16507010      0.21844481      0.15685861      0.14271659
      0.17601876      0.16278913      0.15503383      0.15412144      0.16917585
      0.19107316      0.19426652      0.15822719      0.21479526      0.17647495
      0.14226040      0.18970458      0.22574392      0.20293421      0.21297048
      0.21798862      0.17236921      0.21023332      0.20703996      0.17647495
      0.22072578      0.22802489      0.16689488      0.18103689      0.15822719
      0.21160190      0.20932093      0.15184047      0.21798862      0.14271659
      0.20430279      0.14545376      0.16871965      0.18012451      0.19289794
      0.19107316      0.22026959      0.17100062      0.14955950      0.21388287
      0.18833600      0.22391914      0.22391914      0.16142055      0.21433906
      0.22802489      0.18012451      0.16643868      0.18514264      0.20384660
      0.21981339      0.16187674      0.21114570      0.14180420      0.16780727
      0.22528772      0.15412144      0.16598249      0.22118197      0.19745988
      0.17693115      0.21525145      0.21479526      0.15047189      0.16415771
      0.20567138      0.22483153      0.20886473      0.20977712      0.20840854
      0.16963204      0.22483153      0.15777100      0.21251429      0.18559884
      0.14545376      0.15001570      0.15731480      0.20019705      0.16187674
      0.19244175      0.15184047      0.20430279      0.15503383      0.16735107
      0.19517891      0.22483153      0.17236921      0.17373779      0.20247802
      0.15092808      0.21023332      0.22893727      0.19426652      0.22437533
      0.14317278      0.16324532      0.21023332      0.18787981      0.14180420
      0.16917585      0.20749615      0.22620011      0.23076205      0.15320905
      0.20202182      0.16233294      0.16735107      0.21798862      0.14819092
      0.22255056      0.22391914      0.14955950      0.20612757      0.16643868
      0.16142055      0.19837227      0.16461391      0.16826346      0.17966831
      0.15047189      0.15913958      0.14773473      0.20202182      0.19882846
      0.19517891      0.22437533      0.15412144      0.20384660      0.18879219
      0.19381033      0.15913958      0.20293421      0.21616384      0.18514264
      0.16552630      0.14454137      0.15229667      0.15366525      0.18970458
      0.18924839      0.20202182      0.17419398      0.18696742      0.22848108
      0.21570765      0.16871965      0.19517891      0.22665630      0.19152936
      0.20567138      0.15594622      0.15777100      0.17465018      0.20795235
      0.15001570      0.19426652      0.15229667      0.18879219      0.18423025
      0.20019705      0.15640241      0.16780727      0.21160190      0.15320905
      0.17419398      0.18240548      0.17100062      0.17008824      0.18559884
      0.21707623      0.18286167      0.18879219      0.21205809      0.20521518
      0.21890100      0.16552630      0.15184047      0.14910331      0.15184047
      0.21798862      0.19472272      0.17784354      0.15594622      0.20475899
      0.22346295      0.22756869      0.14819092      0.14454137      0.17875592
      0.15822719      0.22528772      0.21662003      0.14864711      0.20703996
      0.22163817      0.17784354      0.14226040      0.21297048      0.21844481
      0.21251429      0.18787981      0.22802489      0.18240548      0.19609130
      0.14955950      0.21205809      0.18012451      0.16826346      0.20475899
      0.17738734      0.19745988      0.14590995      0.22848108      0.16871965
      0.15503383      0.16963204      0.23030586      0.19016078      0.20567138
      0.16278913      0.17601876      0.22939347      0.20612757      0.16780727
      0.15366525      0.19016078      0.14499756      0.20567138      0.19244175
      0.14727853      0.14682234      0.15640241      0.19517891      0.14590995
      0.22391914      0.20293421      0.18924839      0.21479526      0.20065324
      0.19563511      0.15366525      0.23076205      0.19426652      0.21616384
      0.14271659      0.22300675      0.19426652      0.15184047      0.16735107
      0.18970458      0.21251429      0.15320905      0.14271659      0.15320905
      0.16552630      0.21297048      0.16415771      0.14454137      0.18149309
      0.18559884      0.15229667      0.15047189      0.15959577      0.18423025
      0.19837227      0.17601876      0.14180420      0.14180420      0.20156563
      0.18833600      0.17921212      0.18879219      0.21251429      0.15366525
      0.17236921      0.17647495      0.18833600      0.15366525      0.18012451
      0.19107316      0.14636614      0.17373779      0.17191301      0.21205809
      0.15138428      0.18377406      0.21525145      0.16507010      0.20749615
      0.20749615      0.21935720      0.20202182      0.14590995      0.16552630
      0.21068951      0.14636614      0.21844481      0.22528772      0.20293421
      0.20658376      0.17693115      0.18924839      0.19198555      0.21935720
      0.16142055      0.15366525      0.16187674      0.22026959      0.21844481
      0.19791608      0.16233294      0.15503383      0.21525145      0.19745988
      0.14545376      0.20339040      0.18468645      0.21388287      0.20977712
      0.17784354      0.19107316      0.16963204      0.19016078      0.16826346
      0.17556257      0.19426652      0.19837227      0.14590995      0.20567138
      0.21616384      0.22574392      0.22848108      0.19928466      0.18331786
      0.21981339      0.17373779      0.18605503      0.19928466      0.17373779
      0.20293421      0.21981339      0.23167444      0.18149309      0.17145682
      0.21616384      0.14636614      0.18559884      0.15229667      0.20247802
      0.22437533      0.14590995      0.16917585      0.14910331      0.17373779
      0.15412144      0.17693115      0.20247802      0.16689488      0.17829973
      0.21114570      0.19289794      0.15503383      0.21479526      0.20110943
      0.22437533      0.15685861      0.16552630      0.18559884      0.14819092
      0.22574392      0.15275286      0.15001570      0.21023332      0.18787981
      0.22711250      0.21525145      0.17556257      0.21114570      0.18468645
      0.21525145      0.22711250      0.14910331      0.20567138      0.16278913
      0.22893727      0.19882846      0.14682234      0.17693115      0.18468645
      0.19198555      0.14864711      0.18423025      0.18331786      0.18423025
      0.21297048      0.16735107      0.16826346      0.20247802      0.20612757
      0.16050816      0.17100062      0.17145682      0.22437533      0.19016078
      0.19061697      0.21251429      0.22346295      0.19152936      0.17601876
      0.19289794      0.15184047      0.17465018      0.15092808      0.14773473
      0.15549003      0.16415771      0.20293421      0.19791608      0.19289794
      0.19517891      0.18012451      0.21433906      0.16324532      0.14864711
      0.16233294      0.17236921      0.16187674      0.20384660      0.17145682
      0.18103689      0.21023332      0.15959577      0.17054443      0.15320905
      0.14180420      0.20156563      0.16005197      0.19974085      0.18833600
      0.17921212      0.18286167      0.19791608      0.17419398      0.15640241
      0.19472272      0.19654749      0.21297048      0.17921212      0.17556257
      0.18559884      0.17510637      0.23213063      0.16050816      0.17419398
      0.15001570      0.17328159      0.18605503      0.16689488      0.16643868
      0.19198555      0.21935720      0.17328159      0.21023332      0.18651122
      0.17921212      0.16917585      0.15366525      0.20567138      0.16461391
      0.18970458      0.18331786      0.14271659      0.17647495      0.22346295
      0.16689488      0.19928466      0.21433906      0.17875592      0.16096435
      0.18787981      0.19381033      0.14682234      0.17191301      0.19426652
      0.19289794      0.15138428      0.15275286      0.20384660      0.15549003
      0.15731480      0.18696742      0.14773473      0.19335413      0.19654749
      0.19563511      0.15092808      0.20840854      0.15822719      0.18742361
      0.18696742      0.16187674      0.19016078      0.20339040      0.22346295
      0.14727853      0.14408517      0.17875592      0.20110943      0.20202182
      0.19198555      0.14545376      0.21707623      0.19198555      0.16552630
      0.16278913      0.18605503      0.19107316      0.18696742      0.20612757
      0.18559884      0.22893727      0.17054443      0.20384660      0.23030586
      0.16005197      0.18696742      0.19061697      0.14864711      0.22300675
      0.17784354      0.17738734      0.22620011      0.19198555      0.14362898
      0.18970458      0.19244175      0.15229667      0.17921212      0.17647495
      0.14955950      0.20932093      0.16096435      0.18423025      0.16278913
      0.17875592      0.20521518      0.19517891      0.14362898      0.22665630
      0.17601876      0.20977712      0.18833600      0.18058070      0.20977712
      0.20658376      0.17465018      0.17875592      0.16324532      0.16552630
      0.16689488      0.14545376      0.21525145      0.18194928      0.21890100
      0.17465018      0.20156563      0.16598249      0.20521518      0.20475899
      0.17100062      0.21570765      0.14226040      0.18012451      0.15731480
      0.15959577      0.22528772      0.18970458      0.16963204      0.15275286
      0.16689488      0.20110943      0.15001570      0.15275286      0.20977712
      0.17191301      0.16598249      0.17191301      0.18058070      0.19244175
      0.15913958      0.18970458      0.14727853      0.22711250      0.14362898
      0.16233294      0.23213063      0.21160190      0.19472272      0.16324532
      0.21798862      0.16598249      0.22072578      0.22802489      0.22118197
      0.15685861      0.21525145      0.22346295      0.16826346      0.22026959
      0.16461391      0.15959577      0.14773473      0.21388287      0.19609130
      0.20703996      0.15913958      0.15184047      0.15366525      0.18833600
      0.22939347      0.21160190      0.15366525      0.17875592      0.18194928
      0.14408517      0.19837227      0.20430279      0.18012451      0.14362898
      0.14499756      0.22391914      0.17145682      0.19107316      0.22391914
      0.19107316      0.22802489      0.18468645      0.18833600      0.17236921
      0.18651122      0.21160190      0.15366525      0.22437533      0.18423025
      0.22437533      0.15457764      0.18651122      0.21890100      0.15412144
      0.19974085      0.18605503      0.18696742      0.20019705      0.16415771
      0.21935720      0.15549003      0.18058070      0.20247802      0.15457764
      0.22437533      0.16735107      0.18742361      0.18742361      0.14955950
      0.20475899      0.17875592      0.16689488      0.15822719      0.16050816
      0.21981339      0.15092808      0.17008824      0.19700369      0.17465018
      0.19426652      0.19061697      0.14180420      0.14682234      0.14819092
      0.22574392      0.18468645      0.22665630      0.21068951      0.16507010
      0.20658376      0.22483153      0.14910331      0.15047189      0.16552630
      0.15822719      0.18879219      0.15777100      0.22072578      0.18012451
      0.21251429      0.20156563      0.22255056      0.18696742      0.18012451
      0.20247802      0.15138428      0.17647495      0.15001570      0.18468645
      0.14590995      0.15913958      0.22072578      0.16233294      0.22163817
      0.22574392      0.21753242      0.16278913      0.18194928      0.22848108
      0.16552630      0.15503383      0.23030586      0.21297048      0.17145682
      0.14955950      0.14910331      0.22346295      0.22574392      0.19061697
      0.19609130      0.22665630      0.19198555      0.22026959      0.20065324
      0.18696742      0.22346295      0.21662003      0.16005197      0.21479526
      0.23030586      0.16142055      0.19609130      0.19928466      0.14226040
      0.21981339      0.21890100      0.16871965      0.15685861      0.20293421
      0.15001570      0.19198555      0.20384660      0.16826346      0.17647495
      0.15320905      0.15320905      0.23213063      0.18514264      0.17328159
      0.19472272      0.17510637      0.20795235      0.14408517      0.14408517
      0.21479526      0.14590995      0.16370152      0.18787981      0.16233294
      0.22574392      0.23076205      0.16278913      0.17875592      0.20339040
      0.14819092      0.16598249      0.15959577      0.15913958      0.14590995
      0.16689488      0.14955950      0.19152936      0.17054443      0.22300675
      0.20612757      0.16233294      0.18879219      0.20977712      0.16233294
      0.21160190      0.16233294      0.19381033      0.21753242      0.22665630
      0.21160190      0.21479526      0.22255056      0.20019705      0.21981339
      0.18012451      0.16917585      0.17100062      0.16780727      0.16461391
      0.18514264      0.16735107      0.22665630      0.21616384      0.17054443
      0.22300675      0.20703996      0.17419398      0.16324532      0.21433906
      0.17738734      0.18423025      0.18377406      0.21525145      0.14226040
      0.20749615      0.18696742      0.17100062      0.16324532      0.14590995
      0.22620011      0.15320905      0.20339040      0.19289794      0.19289794
      0.16735107      0.21342668      0.20886473      0.21753242      0.21890100
      0.19152936      0.21251429      0.15138428      0.14454137      0.14226040
      0.14271659      0.18696742      0.18742361      0.23213063      0.21707623
      0.20749615      0.18514264      0.19016078      0.16689488      0.16643868
      0.19882846      0.20247802      0.15913958      0.14454137      0.19791608
      0.16461391      0.19974085      0.14408517      0.20932093      0.19335413
      0.18696742      0.21707623      0.21570765      0.19335413      0.16780727
      0.20567138      0.14636614      0.21707623      0.19837227      0.18103689
      0.21160190      0.17556257      0.19381033      0.14499756      0.19335413
      0.17328159      0.19517891      0.18058070      0.16963204      0.16871965
      0.21798862      0.20247802      0.19928466      0.20156563      0.22756869
      0.14408517      0.15594622      0.21935720      0.20430279      0.18058070
      0.22346295      0.20521518      0.17693115      0.16187674      0.18605503
      0.20202182      0.21433906      0.19061697      0.20658376      0.19244175
      0.19654749      0.21433906      0.16050816      0.14271659      0.20612757
      0.22163817      0.19335413      0.16415771      0.15549003      0.19152936
      0.22072578      0.16005197      0.14317278      0.15594622      0.21890100
      0.21753242      0.16507010      0.16096435      0.21160190      0.14773473
      0.21068951      0.19152936      0.20202182      0.14773473      0.16689488
      0.19745988      0.21479526      0.20795235      0.16233294      0.20019705
      0.19426652      0.15731480      0.22163817      0.16598249      0.21616384
      0.17282540      0.19745988      0.21890100      0.14773473      0.15092808
      0.22346295      0.16917585      0.17008824      0.15184047      0.22574392
      0.21160190      0.15366525      0.16187674      0.14180420      0.18194928
      0.20886473      0.15138428      0.18924839      0.17784354      0.19654749
      0.19381033      0.18149309      0.16780727      0.17556257      0.17236921
      0.22026959      0.22802489      0.22300675      0.20658376      0.14317278
      0.18514264      0.16324532      0.15366525      0.16780727      0.22026959
      0.18103689      0.22756869      0.16917585      0.21525145      0.18423025
      0.22209436      0.15001570      0.21890100      0.22848108      0.17328159
      0.17556257      0.18787981      0.15822719      0.19152936      0.15047189
      0.22665630      0.20065324      0.23167444      0.22802489      0.19381033
      0.16370152      0.19335413      0.16598249      0.19107316      0.14499756
      0.15503383      0.21023332      0.20840854      0.21297048      0.18377406
      0.22939347      0.20019705      0.21251429      0.17966831      0.22163817
      0.19609130      0.17647495      0.20567138      0.16324532      0.15959577
      0.22711250      0.19609130      0.19563511      0.19700369      0.22255056
      0.16278913      0.21753242      0.18103689      0.15092808      0.17693115
      0.16415771      0.17966831      0.20795235      0.23121824      0.19791608
      0.16142055      0.19198555      0.15777100      0.21707623      0.19472272
      0.19974085      0.16278913      0.22939347      0.17236921      0.15549003
      0.20247802      0.17328159      0.22939347      0.18879219      0.22756869
      0.20247802      0.14499756      0.20019705      0.17556257      0.18468645
      0.15229667      0.18058070      0.22391914      0.15640241      0.17647495
      0.16370152      0.15777100      0.18742361      0.15320905      0.21981339
      0.19016078      0.17875592      0.22574392      0.21388287      0.21479526
      0.14955950      0.20293421      0.16871965      0.21525145      0.16826346
      0.14910331      0.22026959      0.14180420      0.23121824      0.21342668
      0.15229667      0.23076205      0.22848108      0.18103689      0.19335413
      0.22756869      0.17601876      0.22255056      0.15184047      0.18012451
      0.20475899      0.22620011      0.18742361      0.15366525      0.20156563
      0.20293421      0.17465018      0.20110943      0.18970458      0.15868338
      0.21570765      0.21297048      0.20749615      0.15777100      0.16461391
      0.19700369      0.20567138      0.21890100      0.15138428      0.14454137
      0.20612757      0.22665630      0.19609130      0.18924839      0.21433906
      0.19974085      0.16735107      0.19974085      0.16963204      0.22711250
      0.22255056      0.20886473      0.20840854      0.14180420      0.17693115
      0.22118197      0.23076205      0.19517891      0.15366525      0.15685861
      0.16826346      0.19472272      0.16142055      0.14545376      0.22848108
      0.16917585      0.21935720      0.18696742      0.15184047      0.20384660
      0.15184047      0.18970458      0.16050816      0.20795235      0.22255056
      0.17784354      0.22665630      0.21433906      0.15640241      0.18605503
      0.22893727      0.16917585      0.22300675      0.18970458      0.17419398
      0.19426652      0.22802489      0.15092808      0.23213063      0.22574392
      0.18696742      0.14271659      0.14773473      0.15092808      0.21297048
      0.17921212      0.16278913      0.17008824      0.19289794      0.17008824
      0.20065324      0.14226040      0.19426652      0.23076205      0.23167444
      0.21935720      0.16826346      0.18058070      0.17966831      0.19152936
      0.22620011      0.20430279      0.22346295      0.17008824      0.17008824
      0.17465018      0.16050816      0.22939347      0.19609130      0.19472272
      0.17373779      0.18377406      0.14226040      0.15412144      0.17875592
      0.21251429      0.17236921      0.22391914      0.17100062      0.17100062
      0.18058070      0.15138428      0.18651122      0.20430279      0.19107316
      0.22574392      0.14180420      0.22391914      0.18605503      0.20749615
      0.15229667      0.14545376      0.22620011      0.20977712      0.16278913
      0.14682234      0.21342668      0.17829973      0.17145682      0.18103689
      0.21023332      0.20475899      0.16050816      0.14271659      0.22893727
      0.15320905      0.15777100      0.15913958      0.21616384      0.19472272
      0.15366525      0.15275286      0.18194928      0.18331786      0.20703996
      0.16233294      0.15822719      0.15184047      0.18240548      0.23030586
      0.21297048      0.17008824      0.22255056      0.18103689      0.17738734
      0.20247802      0.22072578      0.15412144      0.17510637      0.18377406
      0.20019705      0.16324532      0.20339040      0.23121824      0.22163817
      0.20932093      0.14499756      0.22300675      0.21342668      0.15640241
      0.19289794      0.23076205      0.19745988      0.18970458      0.17100062
      0.21753242      0.16598249      0.19381033      0.20156563      0.15366525
      0.16461391      0.22574392      0.20703996      0.16324532      0.22072578
      0.15138428      0.18787981      0.18331786      0.18924839      0.14682234
      0.15092808      0.14910331      0.22483153      0.18468645      0.18149309
      0.15868338      0.22756869      0.19061697      0.18924839      0.22346295
      0.22209436      0.22711250      0.15366525      0.22711250      0.22893727
      0.22346295      0.22209436      0.21753242      0.16370152      0.22255056
      0.22984966      0.19882846      0.15001570      0.17738734      0.20293421
      0.15731480      0.20567138      0.16643868      0.17601876      0.19700369
      0.17829973      0.16050816      0.17236921      0.18331786      0.22848108
      0.15959577      0.14910331      0.19472272      0.15503383      0.21890100
      0.15457764      0.14226040      0.22118197      0.18286167      0.17373779
      0.15777100      0.19563511      0.19882846      0.21068951      0.15685861
      0.15092808      0.22939347      0.17966831      0.15959577      0.22711250
      0.21935720      0.14727853      0.18924839      0.18377406      0.15822719
      0.16917585      0.19609130      0.15594622      0.19061697      0.17191301
      0.23030586      0.18924839      0.17328159      0.17373779      0.19928466
      0.17373779      0.19244175      0.17556257      0.15731480      0.21479526
      0.19152936      0.22209436      0.20156563      0.19563511      0.14362898
      0.15138428      0.17054443      0.19107316      0.18970458      0.19289794
      0.21753242      0.19609130      0.20886473      0.18970458      0.21342668
      0.21753242      0.16826346      0.20703996      0.21935720      0.18696742
      0.21570765      0.23121824      0.16005197      0.19426652      0.15138428
      0.14727853      0.18240548      0.14910331      0.19289794      0.15138428
      0.15913958      0.16005197      0.16324532      0.17784354      0.17282540
      0.19791608      0.16461391      0.17829973      0.19016078      0.17556257
      0.15184047      0.20658376      0.17556257      0.20430279      0.22939347
      0.18012451      0.14545376      0.21160190      0.20612757      0.15047189
      0.20110943      0.23213063      0.20658376      0.21707623      0.14590995
      0.17510637      0.21935720      0.16278913      0.14226040      0.16917585
      0.20703996      0.14362898      0.14271659      0.22072578      0.18103689
      0.23167444      0.18742361      0.20293421      0.18559884      0.19882846
      0.19609130      0.19472272      0.17100062      0.19152936      0.22346295
      0.22346295      0.19609130      0.15366525      0.20110943      0.16324532
      0.22574392      0.16142055      0.18742361      0.17373779      0.14682234
      0.19882846      0.16643868      0.19152936      0.22391914      0.18194928
      0.18286167      0.16871965      0.22665630      0.21890100      0.17373779
      0.21616384      0.15320905      0.21753242      0.22163817      0.22848108
      0.18058070      0.16689488      0.15412144      0.17373779      0.14408517
      0.14408517      0.15549003      0.16735107      0.20795235      0.17191301
      0.20430279      0.14454137      0.19335413      0.19198555      0.22756869
      0.15959577      0.15457764      0.16507010      0.21114570      0.15275286
      0.18559884      0.18149309      0.19152936      0.16598249      0.18696742
      0.19745988      0.18833600      0.14682234      0.18423025      0.17236921
      0.22072578      0.17419398      0.21251429      0.20703996      0.18742361
      0.20795235      0.16552630      0.15047189      0.18331786      0.15092808
      0.19654749      0.21753242      0.18058070      0.19107316      0.20065324
      0.18240548      0.16415771      0.22209436      0.20065324      0.22848108
      0.15275286      0.18970458      0.16507010      0.20658376      0.22893727
      0.14408517      0.16507010      0.16963204      0.20612757      0.18696742
      0.16917585      0.22346295      0.15822719      0.20110943      0.20019705
      0.17465018      0.14819092      0.20703996      0.16507010      0.21479526
      0.18559884      0.19198555      0.15731480      0.20110943      0.17966831
      0.23213063      0.17510637      0.21890100      0.15503383      0.16096435
      0.19016078      0.22437533      0.23030586      0.14226040      0.21525145
      0.23167444      0.16005197      0.18833600      0.18331786      0.17556257
      0.21114570      0.17784354      0.22483153      0.16780727      0.21616384
      0.21023332      0.22528772      0.22391914      0.19381033      0.22118197
      0.16142055      0.14454137      0.21297048      0.15731480      0.15184047
      0.22711250      0.22391914      0.16005197      0.19472272      0.14454137
      0.15457764      0.20430279      0.20384660      0.15777100      0.15868338
      0.20475899      0.20019705      0.22756869      0.17236921      0.18833600
      0.20567138      0.21479526      0.20019705      0.17100062      0.20384660
      0.21525145      0.16963204      0.15640241      0.19609130      0.17419398
      0.16324532      0.17647495      0.14180420      0.21707623      0.14317278
      0.18468645      0.21251429      0.14271659      0.16050816      0.22163817
      0.18286167      0.20521518      0.21388287      0.21114570      0.23030586
      0.21981339      0.15275286      0.15640241      0.18468645      0.16689488
      0.16917585      0.16278913      0.20567138      0.16187674      0.21844481
      0.22756869      0.19107316      0.18742361      0.22802489      0.21388287
      0.18970458      0.23076205      0.17008824      0.18514264      0.19609130
      0.20065324      0.18468645      0.19289794      0.17510637      0.18194928
      0.22665630      0.18012451      0.20384660      0.21433906      0.15731480
      0.21388287      0.20293421      0.18970458      0.22756869      0.20110943
      0.17419398      0.22939347      0.14362898      0.16187674      0.19654749
      0.19016078      0.19928466      0.15868338      0.15184047      0.16142055
      0.17236921      0.21707623      0.14545376      0.15412144      0.21342668
      0.18194928      0.16278913      0.18696742      0.22574392      0.16689488
      0.22848108      0.17966831      0.18423025      0.21981339      0.20156563
      0.20977712      0.20247802      0.16278913      0.14408517      0.18468645
      0.21342668      0.22209436      0.18559884      0.16963204      0.17601876
      0.16552630      0.22255056      0.18514264      0.18924839      0.22893727
      0.18833600      0.16917585      0.18149309      0.22711250      0.22893727
      0.22026959      0.22756869      0.18833600      0.17784354      0.20065324
      0.20019705      0.19791608      0.18970458      0.17145682      0.14180420
      0.21479526      0.20840854      0.15412144      0.14864711      0.17145682
      0.22756869      0.20156563      0.22939347      0.15092808      0.19654749
      0.19837227      0.20977712      0.22711250      0.20065324      0.19426652
      0.14271659      0.17419398      0.20339040      0.15503383      0.22848108
      0.14317278      0.21479526      0.15366525      0.14590995      0.15092808
      0.18970458      0.14819092      0.19016078      0.14545376      0.19654749
      0.21707623      0.17419398      0.15913958      0.16278913      0.15092808
      0.19061697      0.19928466      0.18012451      0.22984966      0.19289794
      0.15640241      0.20795235      0.19609130      0.22665630      0.16598249
      0.16917585      0.17601876      0.17373779      0.18514264      0.17556257
      0.18970458      0.22300675      0.21753242      0.21114570      0.19609130
      0.20384660      0.18286167      0.21707623      0.14454137      0.18514264
      0.15001570      0.18833600      0.20293421      0.20612757      0.17145682
      0.22984966      0.21890100      0.15047189      0.19198555      0.16826346
      0.15138428      0.18331786      0.21662003      0.14955950      0.17328159
      0.16005197      0.18240548      0.14226040      0.21525145      0.22437533
      0.18879219      0.22118197      0.17054443      0.17054443      0.23167444
      0.19882846      0.21753242      0.22939347      0.17465018      0.20521518
      0.19745988      0.22437533      0.21981339      0.16370152      0.16643868
      0.15685861      0.20475899      0.22163817      0.22437533      0.15822719
      0.20156563      0.16826346      0.15320905      0.19289794      0.21114570
      0.15685861      0.20658376      0.21433906      0.17054443      0.19381033
      0.14864711      0.16324532      0.21616384      0.16871965      0.16552630
      0.15229667      0.16643868      0.19791608      0.14864711      0.14180420
      0.22483153      0.22620011      0.17236921      0.18924839      0.18879219
      0.21114570      0.15047189      0.22437533      0.17966831      0.17921212
      0.21297048      0.16826346      0.19974085      0.22118197      0.18696742
      0.15229667      0.22984966      0.14317278      0.17510637      0.20339040
      0.21479526      0.16963204      0.21251429      0.19198555      0.22255056
      0.17191301      0.14545376      0.18833600      0.18559884      0.19381033
      0.20932093      0.16507010      0.14454137      0.16917585      0.15913958
      0.20019705      0.23076205      0.21023332      0.20110943      0.17465018
      0.21433906      0.17510637      0.16507010      0.18787981      0.20339040
      0.20703996      0.19061697      0.17601876      0.18742361      0.22209436
      0.19426652      0.21388287      0.15959577      0.18012451      0.21981339
      0.17829973      0.18696742      0.15001570      0.15001570      0.19472272
      0.23121824      0.21844481      0.15503383      0.20840854      0.16050816
      0.18605503      0.16233294      0.14180420      0.22848108      0.23167444
      0.18559884      0.22620011      0.17008824      0.15959577      0.18240548
      0.16415771      0.14636614      0.19016078      0.14864711      0.16187674
      0.20019705      0.16643868      0.18879219      0.15047189      0.20475899
      0.20795235      0.20977712      0.21753242      0.15001570      0.16142055
      0.15412144      0.18149309      0.17921212      0.16826346      0.20703996
      0.18194928      0.22346295      0.17282540      0.16187674      0.14819092
      0.22209436      0.22802489      0.15092808      0.14317278      0.16963204
      0.17647495      0.17510637      0.15320905      0.14864711      0.16050816
      0.14362898      0.18286167      0.22483153      0.22939347      0.19152936
      0.22893727      0.17100062      0.14955950      0.15366525      0.18103689
      0.14226040      0.16233294      0.14545376      0.17738734      0.16826346
      0.20840854      0.16643868      0.17693115      0.21616384      0.14864711
      0.22255056      0.17647495      0.15731480      0.22163817      0.16233294
      0.16963204      0.22665630      0.14727853      0.17966831      0.14636614
      0.19244175      0.15685861      0.14819092      0.15594622      0.22848108
      0.17373779      0.15777100      0.20886473      0.16142055      0.21160190
      0.18696742      0.19882846      0.20019705      0.18423025      0.20932093
      0.14590995      0.18696742      0.16780727      0.19107316      0.16507010
      0.18742361      0.16735107      0.15412144      0.18240548      0.15229667
      0.20977712      0.14864711      0.14226040      0.14408517      0.16917585
      0.15366525      0.18058070      0.19837227      0.23030586      0.20475899
      0.19700369      0.17510637      0.15184047      0.17556257      0.19517891
      0.16598249      0.15777100      0.14910331      0.22574392      0.21570765
      0.21935720      0.19837227      0.22665630      0.19107316      0.18194928
      0.15594622      0.15777100      0.22893727      0.19335413      0.19061697
      0.14773473      0.18331786      0.18879219      0.16507010      0.22163817
      0.15457764      0.17875592      0.15594622      0.19335413      0.20430279
      0.16507010      0.14408517      0.18696742      0.22437533      0.19426652
      0.19517891      0.18833600      0.20703996      0.21981339      0.21570765
      0.23213063      0.15001570      0.21662003      0.18924839      0.21844481
      0.22209436      0.16643868      0.15503383      0.22437533      0.15184047
      0.16780727      0.23167444      0.19563511      0.17693115      0.18331786
      0.21023332      0.19791608      0.23121824      0.17601876      0.22711250
      0.17921212      0.15868338      0.19654749      0.15777100      0.21616384
      0.15594622      0.20567138      0.14226040      0.15685861      0.21616384
      0.21068951      0.20567138      0.14454137      0.23213063      0.21388287
      0.19745988      0.22391914      0.23167444      0.16324532      0.14499756
      0.22209436      0.17282540      0.21890100      0.15913958      0.21205809
      0.18879219      0.14636614      0.20430279      0.19016078      0.23076205
      0.18605503      0.14819092      0.16461391      0.19016078      0.16461391
      0.22711250      0.18240548      0.20430279      0.20339040      0.15959577
      0.18468645      0.21342668      0.14545376      0.14819092      0.22711250
      0.20521518      0.18377406      0.14590995      0.17921212      0.19061697
      0.18058070      0.20977712      0.14590995      0.17373779      0.22072578
      0.18058070      0.18970458      0.16735107      0.17008824      0.18103689
      0.14271659      0.16871965      0.15640241      0.19517891      0.19335413
      0.17008824      0.22756869      0.23167444      0.22939347      0.15047189
      0.23076205      0.22072578      0.21297048      0.21433906      0.16507010
      0.20567138      0.16507010      0.17784354      0.22026959      0.22893727
      0.23030586      0.19837227      0.21433906      0.17373779      0.18879219
      0.20749615      0.20475899      0.22665630      0.20247802      0.20065324
      0.16278913      0.17419398      0.21890100      0.16826346      0.17556257
      0.15138428      0.14864711      0.22665630      0.15229667      0.17875592
      0.18924839      0.16826346      0.20430279      0.17191301      0.16689488
      0.16598249      0.16461391      0.17008824      0.21570765      0.19107316
      0.16096435      0.22711250      0.15640241      0.16917585      0.20567138
      0.22391914      0.16142055      0.18742361      0.18012451      0.18331786
      0.19517891      0.18331786      0.18559884      0.22391914      0.21251429
      0.15092808      0.15594622      0.21023332      0.15047189      0.21205809
      0.14545376      0.22391914      0.19381033      0.21525145      0.17100062
      0.16735107      0.20521518      0.17921212      0.20977712      0.21342668
      0.17510637      0.17419398      0.15457764      0.23213063      0.17966831
      0.21388287      0.17008824      0.20795235      0.22802489      0.20521518
      0.16963204      0.16917585      0.22255056      0.17145682      0.16507010
      0.22574392      0.18468645      0.18970458      0.21433906      0.16507010
      0.21068951      0.20795235      0.15777100      0.19745988      0.18787981
      0.17373779      0.15913958      0.16689488      0.18514264      0.22209436
      0.18194928      0.23167444      0.23213063      0.18605503      0.19107316
      0.22802489      0.16142055      0.15594622      0.18651122      0.18058070
      0.16826346      0.14271659      0.15549003      0.19974085      0.15229667
      0.16643868      0.21570765      0.14454137      0.15868338      0.15412144
      0.16871965      0.22163817      0.14545376      0.15640241      0.17601876
      0.17966831      0.15685861      0.21114570      0.19974085      0.19198555
      0.15412144      0.20612757      0.16735107      0.14362898      0.15275286
      0.20886473      0.16917585      0.19472272      0.18833600      0.21616384
      0.22209436      0.18149309      0.15685861      0.21023332      0.14362898
      0.17784354      0.17966831      0.16370152      0.20521518      0.21479526
      0.17054443      0.20384660      0.16324532      0.16461391      0.18331786
      0.19517891      0.16370152      0.21297048      0.21023332      0.19928466
      0.19654749      0.22756869      0.17191301      0.18787981      0.23030586
      0.16187674      0.19745988      0.15092808      0.16187674      0.22163817
      0.21570765      0.19198555      0.21433906      0.23213063      0.14408517
      0.22802489      0.16005197      0.21160190      0.19426652      0.21433906
      0.19244175      0.21798862      0.22300675      0.22893727      0.14271659
      0.15913958      0.15366525      0.17100062      0.16917585      0.19517891
      0.15777100      0.19244175      0.18240548      0.19745988      0.18468645
      0.21342668      0.21297048      0.19745988      0.17693115      0.15868338
      0.17556257      0.18377406      0.15640241      0.21160190      0.22346295
      0.14773473      0.14910331      0.20521518      0.16187674      0.16780727
      0.16050816      0.21342668      0.14545376      0.18605503      0.19152936
      0.16598249      0.22848108      0.15503383      0.18742361      0.18605503
      0.19244175      0.18514264      0.17738734      0.16689488      0.18423025
      0.17510637      0.15594622      0.15594622      0.19472272      0.16871965
      0.18651122      0.17100062      0.16963204      0.17100062      0.15275286
      0.17647495      0.18833600      0.17921212      0.23030586      0.22118197
      0.20612757      0.16005197      0.20703996      0.23121824      0.21981339
      0.16735107      0.14727853      0.18787981      0.19609130      0.17647495
      0.21068951      0.14590995      0.20521518      0.19016078      0.18377406
      0.21707623      0.18559884      0.15412144      0.21068951      0.16233294
      0.22391914      0.19061697      0.17328159      0.15275286      0.14910331
      0.20658376      0.16096435      0.18331786      0.16278913      0.18331786
      0.16689488      0.15412144      0.21616384      0.19016078      0.20430279
      0.22939347      0.19244175      0.14317278      0.17875592      0.17465018
      0.16142055      0.20658376      0.18970458      0.14682234      0.21844481
      0.17236921      0.21844481      0.17738734      0.20293421      0.15001570
      0.20977712      0.17054443      0.15685861      0.14362898      0.16826346
      0.20019705      0.22802489      0.22984966      0.22026959      0.16278913
      0.16096435      0.22574392      0.17328159      0.19563511      0.17556257
      0.14955950      0.22939347      0.21570765      0.20475899      0.16370152
      0.15594622      0.16643868      0.18924839      0.21023332      0.17419398
      0.15047189      0.14819092      0.15549003      0.15913958      0.14362898
      0.19244175      0.16689488      0.22391914      0.16963204      0.21479526
      0.14499756      0.22756869      0.14226040      0.17829973      0.18377406
      0.14910331      0.19563511      0.18879219      0.21433906      0.14317278
      0.16871965      0.21798862      0.22391914      0.19016078      0.23167444
      0.23213063      0.21433906      0.18742361      0.20521518      0.21844481
      0.16917585      0.17054443      0.21114570      0.15092808      0.23167444
      0.20110943      0.17236921      0.20293421      0.20612757      0.21844481
      0.15731480      0.14955950      0.14636614      0.14454137      0.18514264
      0.17282540      0.16689488      0.14362898      0.18468645      0.20339040
      0.14955950      0.14955950      0.17921212      0.21890100      0.19381033
      0.20795235      0.21662003      0.23076205      0.18696742      0.17693115
      0.14545376      0.15229667      0.14408517      0.15503383      0.15640241
      0.15913958      0.15229667      0.16780727      0.18194928      0.14408517
      0.20703996      0.14408517      0.23121824      0.16507010      0.19381033
      0.16324532      0.18103689      0.18924839      0.19381033      0.20840854
      0.21114570      0.15457764      0.19974085      0.20384660      0.17510637
      0.21981339      0.14317278      0.18240548      0.16096435      0.21205809
      0.21981339      0.14317278      0.15503383      0.19837227      0.18240548
      0.22802489      0.20840854      0.17601876      0.21662003      0.18742361
      0.14773473      0.19882846      0.19517891      0.17145682      0.20612757
      0.20110943      0.22026959      0.14271659      0.22574392      0.18924839
      0.22939347      0.17829973      0.18103689      0.22072578      0.17738734
      0.19837227      0.20795235      0.19061697      0.22665630      0.14499756
      0.17008824      0.16370152      0.20977712      0.18194928      0.22756869
      0.19107316      0.16187674      0.22848108      0.20247802      0.16278913
      0.19016078      0.18377406      0.20932093      0.17921212      0.18514264
      0.17510637      0.20475899      0.21433906      0.14545376      0.18924839
      0.15777100      0.15503383      0.22756869      0.22574392      0.21662003
      0.22255056      0.21160190      0.19700369      0.15594622      0.19472272
      0.20612757      0.17145682      0.14910331      0.16598249      0.22893727
      0.15138428      0.15001570      0.15503383      0.22620011      0.19335413
      0.17738734      0.20430279      0.20110943      0.18605503      0.21433906
      0.20019705      0.16963204      0.17510637      0.18377406      0.16552630
      0.22802489      0.20293421      0.18559884      0.21479526      0.16142055
      0.20156563      0.16552630      0.16005197      0.17921212      0.18787981
      0.16598249      0.21570765      0.20886473      0.18651122      0.16598249
      0.15822719      0.18559884      0.16917585      0.20658376      0.22346295
      0.17236921      0.22893727      0.15457764      0.15047189      0.17236921
      0.19791608      0.14317278      0.20932093      0.18879219      0.19426652
      0.18286167      0.18742361      0.23213063      0.17556257      0.21388287
      0.17145682      0.19016078      0.16370152      0.21935720      0.22163817
      0.18468645      0.19837227      0.18377406      0.20977712      0.17008824
      0.18103689      0.20293421      0.17784354      0.19244175      0.17829973
      0.15685861      0.22939347      0.17875592      0.14910331      0.21297048
      0.19335413      0.17921212      0.16917585      0.16005197      0.17966831
      0.21297048      0.17373779      0.21388287      0.17191301      0.16963204
      0.18377406      0.19700369      0.21251429      0.14499756      0.14864711
      0.18514264      0.23030586      0.17419398      0.18194928      0.15457764
      0.16050816      0.22072578      0.18696742      0.20521518      0.15320905
      0.18879219      0.20567138      0.22620011      0.22391914      0.19837227
      0.18468645      0.22255056      0.17510637      0.17738734      0.19517891
      0.17328159      0.16278913      0.18286167      0.22346295      0.20339040
      0.22574392      0.20612757      0.14408517      0.16552630      0.19609130
      0.16735107      0.15229667      0.20886473      0.15001570      0.20567138
      0.18833600      0.15412144      0.17601876      0.21388287      0.16552630
      0.15685861      0.17465018      0.18468645      0.18970458      0.18879219
      0.17875592      0.19882846      0.18787981      0.19563511      0.17328159
      0.21479526      0.20977712      0.21981339      0.18468645      0.23213063
      0.22574392      0.16552630      0.16461391      0.16370152      0.14955950
      0.18514264      0.20065324      0.21023332      0.22483153      0.21890100
      0.16507010      0.20384660      0.20612757      0.19244175      0.22939347
      0.16142055      0.18012451      0.22665630      0.19882846      0.15047189
      0.19745988      0.22802489      0.18559884      0.15731480      0.15777100
      0.20977712      0.22026959      0.21798862      0.21935720      0.16005197
      0.17966831      0.19198555      0.21433906      0.18194928      0.16233294
      0.21753242      0.19974085      0.18103689      0.22209436      0.17784354
      0.16278913      0.16552630      0.21114570      0.23121824      0.18924839
      0.21935720      0.16598249      0.15138428      0.23030586      0.15503383
      0.16735107      0.19517891      0.21068951      0.22893727      0.17601876
      0.20339040      0.20247802      0.20932093      0.22802489      0.16187674
      0.19974085      0.20247802      0.14819092      0.17693115      0.15320905
      0.16826346      0.22437533      0.19745988      0.21890100      0.20430279
      0.19198555      0.18879219      0.15138428      0.21023332      0.19472272
      0.15092808      0.15731480      0.17647495      0.21662003      0.14590995
      0.16552630      0.17054443      0.19107316      0.18194928      0.16598249
      0.14408517      0.16461391      0.15412144      0.18377406      0.16963204
      0.18012451      0.15959577      0.14362898      0.15549003      0.20521518
      0.20658376      0.16005197      0.15412144      0.16461391      0.14682234
      0.22984966      0.21981339      0.19563511      0.15822719      0.19563511
      0.19928466      0.14499756      0.14499756      0.16598249      0.17829973
      0.15685861      0.20795235      0.19974085      0.15868338      0.14682234
      0.22391914      0.14362898      0.14955950      0.22255056      0.15275286
      0.17008824      0.20749615      0.21023332      0.22346295      0.17145682
      0.16187674      0.16963204      0.20840854      0.15640241      0.18833600
      0.15777100      0.21525145      0.16096435      0.20658376      0.20293421
      0.22939347      0.17738734      0.20521518      0.19928466      0.18787981
      0.19745988      0.19198555      0.19335413      0.18696742      0.19974085
      0.23076205      0.22026959      0.18879219      0.22118197      0.19563511
      0.16278913      0.20521518      0.14955950      0.18286167      0.14180420
      0.21205809      0.18970458      0.21798862      0.18696742      0.20795235
      0.20977712      0.20247802      0.14727853      0.17008824      0.22072578
      0.21342668      0.19198555      0.18696742      0.18331786      0.19198555
      0.17738734      0.20019705      0.16780727      0.16507010      0.22802489
      0.20840854      0.18559884      0.15138428      0.20886473      0.16142055
      0.18468645      0.18012451      0.14362898      0.18833600      0.20065324
      0.15047189      0.19517891      0.22072578      0.21525145      0.22939347
      0.14636614      0.17145682      0.21525145      0.20749615      0.14864711
      0.17100062      0.20658376      0.19016078      0.17738734      0.17966831
      0.20384660      0.20795235      0.21023332      0.18924839      0.18879219
      0.22255056      0.21844481      0.20110943      0.20795235      0.19928466
      0.16005197      0.18696742      0.23121824      0.18605503      0.21616384
      0.22209436      0.17282540      0.21114570      0.18058070      0.20475899
      0.18331786      0.15366525      0.19563511      0.16917585      0.16187674
      0.15320905      0.19381033      0.21068951      0.17100062      0.15868338
      0.18194928      0.15320905      0.19426652      0.19928466      0.20202182
      0.22255056      0.19974085      0.20840854      0.17921212      0.14545376
      0.17510637      0.22255056      0.20475899      0.15047189      0.15868338
      0.16871965      0.16552630      0.18696742      0.21114570      0.20977712
      0.20521518      0.15184047      0.14773473      0.19882846      0.16598249
      0.18696742      0.16096435      0.17556257      0.22026959      0.14819092
      0.15229667      0.14226040      0.21662003      0.16552630      0.20977712
      0.20932093      0.17921212      0.20156563      0.17921212      0.16552630
      0.20612757      0.21981339      0.16415771      0.19016078      0.21388287
      0.18514264      0.18514264      0.19654749      0.17647495      0.17693115
      0.15959577      0.15320905      0.19609130      0.18103689      0.19289794
      0.16050816      0.15229667      0.16552630      0.21570765      0.21616384
      0.21433906      0.22209436      0.16324532      0.14408517      0.22437533
      0.19152936      0.19517891      0.22209436      0.15640241      0.22391914
      0.17465018      0.17054443      0.21160190      0.17054443      0.21433906
      0.20886473      0.18970458      0.19654749      0.19654749      0.21479526
      0.14819092      0.20703996      0.18103689      0.21114570      0.19745988
      0.18331786      0.17784354      0.15229667      0.15412144      0.20475899
      0.16005197      0.15001570      0.17921212      0.18696742      0.15685861
      0.21753242      0.21662003      0.17054443      0.17738734      0.16187674
      0.18286167      0.21616384      0.17738734      0.16187674      0.16142055
      0.20703996      0.17647495      0.18879219      0.15594622      0.16780727
      0.16780727      0.17236921      0.22802489      0.22255056      0.16871965
      0.19061697      0.14362898      0.23167444      0.19016078      0.21890100
      0.21616384      0.20019705      0.22756869      0.15594622      0.15868338
      0.19517891      0.20567138      0.18970458      0.19198555      0.21935720
      0.15320905      0.14727853      0.22620011      0.21388287      0.18514264
      0.18879219      0.21433906      0.20247802      0.14590995      0.15184047
      0.20019705      0.17236921      0.18012451      0.14271659      0.19609130
      0.15822719      0.19928466      0.16507010      0.20202182      0.15503383
      0.15777100      0.20932093      0.22574392      0.14317278      0.15412144
      0.22300675      0.22939347      0.16461391      0.16050816      0.19654749
      0.15184047      0.17100062      0.18651122      0.19563511      0.17008824
      0.15913958      0.15959577      0.21798862      0.20703996      0.19107316
      0.15913958      0.15275286      0.16780727      0.16552630      0.22802489
      0.16643868      0.18058070      0.20475899      0.21753242      0.19472272
      0.19061697      0.17647495      0.20384660      0.16187674      0.22437533
      0.16461391      0.21114570      0.22939347      0.22939347      0.19837227
      0.17829973      0.16780727      0.17100062      0.22665630      0.18468645
      0.16871965      0.16187674      0.22209436      0.22255056      0.14819092
      0.22346295      0.15275286      0.14271659      0.15685861      0.21890100
      0.21433906      0.15822719      0.20156563      0.17601876      0.20795235
      0.19426652      0.18377406      0.22072578      0.19016078      0.14955950
      0.21662003      0.18696742      0.14226040      0.17145682      0.18286167
      0.17693115      0.15320905      0.16598249      0.22072578      0.16643868
      0.15731480      0.20293421      0.16461391      0.20293421      0.18286167
      0.21023332      0.18970458      0.16598249      0.17784354      0.14636614
      0.17282540      0.17647495      0.22391914      0.18058070      0.15731480
      0.14499756      0.21023332      0.20202182      0.21068951      0.20293421
      0.14545376      0.14819092      0.19928466      0.17556257      0.18742361
      0.21844481      0.22848108      0.16415771      0.15913958      0.17465018
      0.17236921      0.16917585      0.15320905      0.20612757      0.19198555
      0.15184047      0.22209436      0.22346295      0.18331786      0.17738734
      0.14317278      0.15047189      0.19198555      0.17556257      0.22848108
      0.18970458      0.18605503      0.14910331      0.22528772      0.22620011
      0.19107316      0.17100062      0.16917585      0.20703996      0.20065324
      0.19517891      0.16187674      0.16871965      0.15092808      0.15503383
      0.15457764      0.19061697      0.20384660      0.17601876      0.14727853
      0.17601876      0.22984966      0.17829973      0.20612757      0.20521518
      0.17054443      0.17875592      0.19472272      0.15184047      0.17419398
      0.20703996      0.18468645      0.17008824      0.16461391      0.20795235
      0.23213063      0.17100062      0.14545376      0.22483153      0.22026959
      0.21479526      0.23213063      0.20339040      0.15320905      0.22574392
      0.18423025      0.15777100      0.21342668      0.21205809      0.20932093
      0.22620011      0.17191301      0.22255056      0.14545376      0.18286167
      0.19563511      0.19244175      0.20384660      0.15503383      0.20567138
      0.14317278      0.17510637      0.21114570      0.19517891      0.22255056
      0.19974085      0.21570765      0.18787981      0.20430279      0.18423025
      0.21844481      0.19745988      0.21890100      0.23121824      0.20110943
      0.17601876      0.14271659      0.16050816      0.21342668      0.14226040
      0.19609130      0.19016078      0.18149309      0.22209436      0.21388287
      0.19198555      0.22848108      0.21890100      0.20293421      0.16233294
      0.16324532      0.19745988      0.19244175      0.19381033      0.19745988
      0.19061697      0.19381033      0.21479526      0.17465018      0.22026959
      0.19791608      0.15959577      0.19061697      0.22346295      0.21023332
      0.16461391      0.18012451      0.22756869      0.20019705      0.22574392
      0.21707623      0.22391914      0.17510637      0.21388287      0.18012451
      0.16826346      0.16871965      0.20749615      0.15685861      0.22346295
      0.16278913      0.22346295      0.23076205      0.21662003      0.14317278
      0.22574392      0.17693115      0.14226040      0.20612757      0.22893727
      0.18423025      0.19745988      0.14499756      0.20293421      0.15047189
      0.21068951      0.20612757      0.19882846      0.22756869      0.18787981
      0.20384660      0.17693115      0.14727853      0.18012451      0.15184047
      0.17647495      0.21616384      0.21662003      0.20202182      0.16871965
      0.16689488      0.21890100      0.19654749      0.17054443      0.23076205
      0.17054443      0.15640241      0.20521518      0.17419398      0.14864711
      0.19289794      0.20293421      0.23076205      0.15503383      0.19974085
      0.16963204      0.16461391      0.20703996      0.17966831      0.22163817
      0.18605503      0.15777100      0.15959577      0.18559884      0.17784354
      0.22711250      0.22711250      0.14180420      0.16963204      0.19517891
      0.23213063      0.14180420      0.20795235      0.20886473      0.17693115
      0.20019705      0.22483153      0.20293421      0.19472272      0.18559884
      0.19244175      0.19745988      0.18742361      0.14636614      0.16598249
      0.14590995      0.17966831      0.21479526      0.20886473      0.14317278
      0.20886473      0.16415771      0.16552630      0.17465018      0.20977712
      0.21342668      0.22939347      0.21844481      0.20430279      0.19837227
      0.21935720      0.14864711      0.17784354      0.17373779      0.21068951
      0.20840854      0.14864711      0.20065324      0.19198555      0.20658376
      0.19882846      0.21251429      0.16050816      0.20475899      0.21068951
      0.22711250      0.20886473      0.17601876      0.18833600      0.19289794
      0.16324532      0.15549003      0.21433906      0.22163817      0.15594622
      0.16415771      0.17601876      0.19381033      0.16005197      0.20247802
      0.14819092      0.20247802      0.15412144      0.14362898      0.18377406
      0.21798862      0.22574392      0.19745988      0.14499756      0.20795235
      0.14819092      0.17647495      0.16963204      0.17966831      0.17738734
      0.22984966      0.20521518      0.22483153      0.23213063      0.14317278
      0.20019705      0.18605503      0.17875592      0.19563511      0.17191301
      0.14408517      0.17054443      0.22528772      0.17784354      0.16050816
      0.14317278      0.20019705      0.18194928      0.22026959      0.19563511
      0.17647495      0.18651122      0.15001570      0.21114570      0.14545376
      0.18377406      0.18559884      0.16871965      0.14271659      0.18879219
      0.21981339      0.15685861      0.14226040      0.21798862      0.22711250
      0.22483153      0.21297048      0.18194928      0.19563511      0.17647495
      0.19609130      0.18787981      0.22528772      0.20339040      0.14499756
      0.19335413      0.17784354      0.18240548      0.15412144      0.21479526
      0.14180420      0.19745988      0.17601876      0.15184047      0.16689488
      0.17921212      0.22255056      0.14636614      0.16096435      0.19882846
      0.16415771      0.22984966      0.22118197      0.14682234      0.16050816
      0.21890100      0.18194928      0.20110943      0.17829973      0.17921212
      0.16233294      0.18377406      0.17465018      0.20202182      0.16735107
      0.14682234      0.21297048      0.15594622      0.17008824      0.15229667
      0.14271659      0.15822719      0.18879219      0.19791608      0.22984966
      0.18103689      0.16963204      0.16689488      0.18970458      0.21068951
      0.19472272      0.21753242      0.21798862      0.22483153      0.17465018
      0.15275286      0.16552630      0.17008824      0.20932093      0.20612757
      0.21753242      0.18787981      0.22939347      0.14317278      0.20521518
      0.17465018      0.17784354      0.19517891      0.14180420      0.18924839
      0.19974085      0.17966831      0.22665630      0.19563511      0.18970458
      0.14636614      0.18468645      0.20019705      0.17419398      0.20612757
      0.16005197      0.21570765      0.14271659      0.21114570      0.16278913
      0.18605503      0.20065324      0.21753242      0.20065324      0.14271659
      0.20475899      0.20521518      0.16050816      0.16096435      0.14864711
      0.22711250      0.18377406      0.16142055      0.20475899      0.20612757
      0.20840854      0.14864711      0.14819092      0.17328159      0.16050816
      0.22163817      0.15868338      0.22574392      0.18423025      0.21023332
      0.19928466      0.22848108      0.22209436      0.15777100      0.18286167
      0.23213063      0.20110943      0.17966831      0.21342668      0.14955950
      0.17373779      0.14590995      0.15777100      0.19335413      0.23121824
      0.19700369      0.19700369      0.15777100      0.19335413      0.19700369
      0.15868338      0.14773473      0.18696742      0.15959577      0.14773473
      0.19563511      0.15001570      0.17328159      0.21525145      0.19517891
      0.17373779      0.23030586      0.14590995      0.14773473      0.22893727
      0.21662003      0.20156563      0.14226040      0.14910331      0.21023332
      0.21023332      0.22620011      0.21525145      0.16963204      0.14819092
      0.22255056      0.19107316      0.14454137      0.18879219      0.15822719
      0.21662003      0.16324532      0.20749615      0.20886473      0.19381033
      0.15047189      0.21981339      0.16552630      0.19609130      0.15731480
      0.18879219      0.22711250      0.16142055      0.19426652      0.20749615
      0.18331786      0.18924839      0.20932093      0.22848108      0.20247802
      0.18514264      0.15047189      0.17647495      0.21890100      0.17465018
      0.16278913      0.16005197      0.15138428      0.22163817      0.20567138
      0.19381033      0.15549003      0.18377406      0.15275286      0.15092808
      0.19928466      0.21844481      0.14271659      0.15184047      0.16871965
      0.15594622      0.20156563      0.15047189      0.19974085      0.14955950
      0.18833600      0.19061697      0.22209436      0.22026959      0.14864711
      0.20019705      0.22665630      0.20658376      0.20612757      0.21205809
      0.16415771      0.21798862      0.15959577      0.16780727      0.20795235
      0.19745988      0.14910331      0.19837227      0.17419398      0.20658376
      0.19061697      0.20886473      0.22939347      0.23167444      0.15868338
      0.19837227      0.17282540      0.16324532      0.22756869      0.22528772
      0.20247802      0.20430279      0.17693115      0.19244175      0.16689488
      0.16963204      0.21570765      0.20293421      0.19152936      0.18149309
      0.15549003      0.19016078      0.21479526      0.23121824      0.19791608
      0.22026959      0.15457764      0.14408517      0.18742361      0.19563511
      0.20612757      0.23121824      0.22528772      0.16917585      0.18286167
      0.20475899      0.19654749      0.18696742      0.19654749      0.14454137
      0.14545376      0.17008824      0.14636614      0.19335413      0.18240548
      0.20065324      0.15001570      0.16142055      0.15412144      0.22209436
      0.16005197      0.20703996      0.16278913      0.15138428      0.15594622
      0.23121824      0.19198555      0.15594622      0.16096435      0.15868338
      0.17282540      0.18605503      0.16871965      0.22574392      0.17419398
      0.19244175      0.15640241      0.21251429      0.21479526      0.15366525
      0.22620011      0.15503383      0.16552630      0.23213063      0.21890100
      0.21068951      0.22483153      0.23076205      0.17829973      0.20886473
      0.14682234      0.19107316      0.15640241      0.18924839      0.19974085
      0.18696742      0.22711250      0.21160190      0.19882846      0.17829973
      0.22984966      0.15001570      0.14271659      0.14864711      0.16507010
      0.15822719      0.17191301      0.18103689      0.16963204      0.14727853
      0.17647495      0.18103689      0.18559884      0.14955950      0.22300675
      0.18149309      0.15503383      0.18970458      0.21935720      0.15184047
      0.19974085      0.15047189      0.22026959      0.22984966      0.21707623
      0.17328159      0.15092808      0.16005197      0.19837227      0.15184047
      0.19107316      0.18331786      0.22711250      0.14819092      0.15731480
      0.20202182      0.19609130      0.21023332      0.18468645      0.16507010
      0.21479526      0.23213063      0.14819092      0.23076205      0.18970458
      0.22437533      0.17282540      0.19152936      0.15275286      0.15229667
      0.18012451      0.18514264      0.19107316      0.19289794      0.18423025
      0.15868338      0.21935720      0.16278913      0.15685861      0.21935720
      0.18058070      0.16187674      0.22939347      0.15549003      0.16871965
      0.14454137      0.18468645      0.19745988      0.15184047      0.16461391
      0.19244175      0.22802489      0.15685861      0.15457764      0.15640241
      0.19609130      0.21935720      0.15320905      0.17008824      0.16005197
      0.22802489      0.15184047      0.20475899      0.23167444      0.18879219
      0.16461391      0.20612757      0.17738734      0.23167444      0.22756869
      0.16826346      0.20521518      0.19244175      0.15549003      0.16278913
      0.22437533      0.22209436      0.16096435      0.20110943      0.21342668
      0.15001570      0.20612757      0.22665630      0.23121824      0.19837227
      0.21479526      0.16050816      0.14362898      0.16233294      0.22984966
      0.20384660      0.22118197      0.17693115      0.21753242      0.20886473
      0.17373779      0.14499756      0.20247802      0.18194928      0.18879219
      0.14499756      0.16142055      0.23167444      0.17784354      0.20247802
      0.19791608      0.19152936      0.15275286      0.17921212      0.16050816
      0.17419398      0.22391914      0.14819092      0.22163817      0.15138428
      0.14864711      0.21388287      0.14727853      0.19198555      0.17510637
      0.21205809      0.15959577      0.22620011      0.14590995      0.16096435
      0.15184047      0.22665630      0.15366525      0.20795235      0.20840854
      0.15549003      0.19335413      0.17054443      0.19152936      0.21662003
      0.15047189      0.19335413      0.14819092      0.22802489      0.14545376
      0.14955950      0.19837227      0.20658376      0.19426652      0.18651122
      0.14727853      0.15959577      0.18879219      0.15731480      0.20749615
      0.17008824      0.20932093      0.20339040      0.23030586      0.20339040
      0.20247802      0.14819092      0.16507010      0.22072578      0.15320905
      0.15503383      0.17236921      0.17510637      0.18879219      0.15092808
      0.16917585      0.17693115      0.14226040      0.19244175      0.14955950
      0.20110943      0.22118197      0.18103689      0.21205809      0.22209436
      0.16598249      0.16871965      0.16735107      0.20339040      0.22528772
      0.18468645      0.16461391      0.23167444      0.21890100      0.20202182
      0.21844481      0.18696742      0.19791608      0.17829973      0.21616384
      0.16735107      0.20749615      0.22072578      0.22437533      0.16689488
      0.22255056      0.18149309      0.17328159      0.20840854      0.19244175
      0.14408517      0.17236921      0.22939347      0.14545376      0.22756869
      0.18423025      0.18468645      0.19745988      0.21160190      0.20019705
      0.14317278      0.19974085      0.20521518      0.17465018      0.16689488
      0.17008824      0.14955950      0.19974085      0.18240548      0.19289794
      0.16187674      0.18468645      0.20156563      0.18012451      0.22300675
      0.22072578      0.15731480      0.23121824      0.14545376      0.16552630
      0.16826346      0.19244175      0.22711250      0.16507010      0.18879219
      0.21570765      0.22939347      0.16871965      0.15640241      0.18559884
      0.15184047      0.19745988      0.16233294      0.17693115      0.14545376
      0.18058070      0.19289794      0.17784354      0.17921212      0.20703996
      0.21205809      0.18058070      0.15184047      0.17328159      0.19700369
      0.18924839      0.18970458      0.21844481      0.21935720      0.20612757
      0.20019705      0.21388287      0.20795235      0.20475899      0.21388287
      0.19882846      0.14499756      0.16461391      0.15913958      0.18103689
      0.21433906      0.22756869      0.17145682      0.20384660      0.21525145
      0.15457764      0.14545376      0.20612757      0.17829973      0.19472272
      0.19745988      0.15868338      0.20156563      0.16689488      0.20932093
      0.22802489      0.15594622      0.20110943      0.22939347      0.14910331
      0.20475899      0.21616384      0.19061697      0.19244175      0.19244175
      0.18194928      0.19152936      0.18696742      0.16415771      0.16324532
      0.18194928      0.20156563      0.19472272      0.17100062      0.15138428
      0.16780727      0.18833600      0.15457764      0.15822719      0.22255056
      0.19974085      0.18423025      0.22574392      0.14226040      0.17556257
      0.16552630      0.18103689      0.19152936      0.18103689      0.16142055
      0.15092808      0.20932093      0.14864711      0.18833600      0.21342668
      0.16598249      0.16643868      0.20977712      0.18058070      0.15092808];
  s = s';
  s = s(:);
end
