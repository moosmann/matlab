% P05 reconstruction pipeline: preprocessing, filtering, phase retrieval,
% tomographic reconstruction, simple segmentation, etc.
%
% USAGE
% Set parameters in PARAMETERS / SETTINGS section below and run script.
%
% HOW TO RUN THE SCRIPT:
% - Editor windows: press 'F5' when focus is in the Editor window
% - Editor tab: click 'Run' in the
% - Command Window: type 'p05_reco' and Enter
%
% HOW TO AUTOMATICALLY LOOP OVER RECONSTRUCTIONS:
% To loop over different data or parameters sets see 'p05_reco_loop_template'.
%
% For some additional information see 'p05_notes'.
%
% Known errors:
% GPU CUDA error: device busy
% Do not call 'gpuDevice()' when ASTRA toolbox is used with multi-GPU
% support, or MATLAB terminates abnormally with a segmentation violation.
%
% Written by Julian Moosmann. First version: 2016-09-28. Last modifcation:
% 2017-11-26

%clear all
close all hidden % close all open windows
%dbstop if error

%% PARAMETERS / SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scan_path = ...
    '/asap3/petra3/gpfs/p05/2017/data/11003288/raw/syn154_58L_Mg_12_cmos_test';
    '/asap3/petra3/gpfs/p05/2017/data/11003773/raw/syn132_cor_mg5gd434s_mg10gd408s_12';
    '/asap3/petra3/gpfs/p05/2017/data/11003288/raw/syn165_58L_Mg10Gd_12w_000'; %reconstruction is crap
    '/asap3/petra3/gpfs/p05/2017/data/11003288/raw/syn166_104R_Mg10Gd_4w_000'; %reconstruction is crap
    '/asap3/petra3/gpfs/p05/2017/data/11003288/raw/syn165_58L_Mg10Gd_12w_003';
'/asap3/petra3/gpfs/p05/2017/data/11003288/raw/syn151_58L_Mg_12_000';
'/asap3/petra3/gpfs/p05/2017/data/11003288/raw/syn154_58L_Mg_12_cmos_test';
'/asap3/petra3/gpfs/p05/2017/data/11003288/raw/syn146_58L_Mg_12_cmos_test';
'/asap3/petra3/gpfs/p05/2017/data/11003288/raw/syn145_58L_Mg_12_cmos_test';
'/asap3/petra3/gpfs/p05/2017/data/11003773/raw/syn133_cor_mg5gd1p_12';
'/asap3/petra3/gpfs/p05/2017/data/11003773/raw/syn131_cor_mg5gd428s_mg10gd407s_mg5gd8p_12';
'/asap3/petra3/gpfs/p05/2017/data/11003773/raw/syn130_cor_mg5gd416s_mg10gd410s_mg5gd7p_12';
'/asap3/petra3/gpfs/p05/2017/data/11003773/raw/syn129_cor_mg5gd413s_mg10gd409s_mg5gd3p_12';
'/asap3/petra3/gpfs/p05/2017/data/11003773/raw/syn128_cor_mg5gd432s_mg10gd418s_pmg10p_12';
'/asap3/petra3/gpfs/p05/2017/data/11003773/raw/syn127_cor_mg5gd430s_mg10gd413s_pmg9p_12';
'/asap3/petra3/gpfs/p05/2017/data/11003656/raw/szeb_19_01';
'/asap3/petra3/gpfs/p05/2017/commissioning/c20171115_000_kit_test/raw/test23_gain60_nb_occddist600mm_400proj';
'/asap3/petra3/gpfs/p05/2017/commissioning/c20171115_000_kit_test/raw/test16_gain60_nb_occddist260mm';
'/asap3/petra3/gpfs/p05/2017/commissioning/c20171115_000_kit_test/raw/test17_gain60_wwb_occddist260mm';
'/asap3/petra3/gpfs/p05/2017/commissioning/c20171115_000_kit_test/raw/test20_gain60_wwbnf_occddist260mm';
'/asap3/petra3/gpfs/p05/2017/commissioning/c20171115_000_kit_test/raw/test21_gain60_wwbnf_occddist260mm_600ref';
'/asap3/petra3/gpfs/p05/2017/commissioning/c20171115_000_kit_test/raw/test17_gain60_wwb_occddist600mm';
'/asap3/petra3/gpfs/p05/2017/commissioning/c20171115_000_kit_test/raw/test18_gain60_nb_occddist600mm';
'/asap3/petra3/gpfs/p05/2017/commissioning/c20171115_000_kit_test/raw/test19_gain60_wwbnf_occddist600mm';
'/asap3/petra3/gpfs/p05/2017/commissioning/c20171115_000_kit_test/raw/test22_gain60_wwbnf_occddist600mm_400proj';
'/asap3/petra3/gpfs/p05/2017/commissioning/c20171115_000_kit_test/raw/test13_gain61_wwb_occddist600mm';
'/asap3/petra3/gpfs/p05/2017/commissioning/c20171115_000_kit_test/raw/test13_gain61_wwb_occddist600mm';
'/asap3/petra3/gpfs/p05/2017/data/11003656/raw/szeb_01';
'/asap3/petra3/gpfs/p05/2017/commissioning/c20171115_000_kit_test/raw/test11_gain60';
'/asap3/petra3/gpfs/p05/2017/commissioning/c20171115_000_kit_test/raw/test10';
'/asap3/petra3/gpfs/p05/2017/commissioning/c20171115_000_kit_test/raw/test9';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn99_43R';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn98_77L';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn96_82L';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn97_74L';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn98_77L';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn99_43R';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn77_cor_P2_2';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn76_cor_P2_1';
%'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn74_cor_P1_6';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn69_cor_P1_I';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn72_cor_P1_4';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn73_cor_P1_5';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn75_cor_P1_7';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn71_cor_P1_3';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn70_cor_P1_2';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn64_CPDexplant';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn65_94R_Mg5Gd_8w';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn67_cor_Punkt1_previous_is_Punkt2';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn68_cor_Punkt3';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn66_cor_Punkt1';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn64_CPDexplant';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn60_cor';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn52_95R_Mg5Gd_8w';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn51_87R_Mg5Gd_4w';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn50_99L_Mg5Gd_4w';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn49_80L_Mg5Gd_8w';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn48_89L_Mg5Gd_4w';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn47_100AL_Mg5Gd_4w';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn46_88R_Mg5Gd_4w';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn45_101BL_Mg5Gd_4w';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn44_66L_Mg5Gd_12w';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn43_38L_PEEK_8w';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn42_38L_PEEK_8w/';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn41_63L_Mg5Gd_12w';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn40_69L_Mg10Gd_12w';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn39_75L_Mg5Gd_8w';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn38_73R_Mg10Gd_8w';
% not all data copied '/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn34_79R_Mg10Gd_8w';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn37_69L_Mg10Gd_12w';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn34_79R_Mg10Gd_8w';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn33_80R_Mg10Gd_8w';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn32_99R_Mg10Gd_4w';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn31_97R_Mg10Gd_4w';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn36_63L_Mg5Gd_12w';
'/asap3/petra3/gpfs/p05/2017/data/11003440/raw/syn35_77R_Mg10Gd_8w';
'/asap3/petra3/gpfs/p05/2017/data/11003700/raw/mpimm_12_a';
'/asap3/petra3/gpfs/p05/2017/data/11003700/raw/mpimm_07_a';
'/asap3/petra3/gpfs/p05/2017/data/11003700/raw/mpimm_04_a';
'/asap3/petra3/gpfs/p05/2017/data/11003435/raw/ony_24'; % fly scan
'/asap3/petra3/gpfs/p05/2017/data/11003435/raw/ony_42';
'/asap3/petra3/gpfs/p05/2017/data/11003950/raw/syn13_55L_Mg10Gd_12w_load_00';
'/asap3/petra3/gpfs/p05/2017/data/11002845/raw/ste_02_l1_bb';
'/asap3/petra3/gpfs/p05/2017/data/11002839/raw/ehh_2017_019_b';
'/asap3/petra3/gpfs/p05/2017/data/11003063/raw/hnee_01_hw_hk776_bn161514';
'/asap3/petra3/gpfs/p05/2017/data/11002839/raw/ehh_2017_015_a';
'/asap3/petra3/gpfs/p05/2017/data/11003950/raw/syn01_48L_PEEK_12w_b';
read_flatcor = 0; % read flatfield-corrected images from disc, skips preprocessing
read_flatcor_path = ''; % subfolder of 'flat_corrected' containing projections
% PREPROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_roi = []; % [y0 y1] vertical roi.  skips first raw_roi(1)-1 lines, reads until raw_roi(2)
raw_bin = 2; % projection binning factor: 1, 2, or 4
bin_before_filtering = 1; % Binning before pixel filtering; much faster but worse filtering
excentric_rot_axis = 0; % off-centered rotation axis increasing FOV. -1: left, 0: centeerd, 1: right. influences rot_corr_area1
crop_at_rot_axis = 0; % for recos of scans with excentric rotation axis but WITHOUT projection stitching
stitch_projections = 0; % for 2 pi scans: stitch projection at rotation axis position
stitch_method = 'sine';'linear'; 'step'; %  ! adjust correlation area if necessary !
% 'step' : no interpolation, use step function
% 'linear' : linear interpolation of overlap region
% 'sine' : sinusoidal interpolation of overlap region
proj_range = 2; % range of projections to be used (from all found). if empty or 1: all, if scalar: stride
ref_range = 2;[]; % range of flat fields to be used (from all found). start:incr:end. if empty or 1: all. if scalar: stride
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
ring_filter = 1; % ring artifact filter
ring_filter_method = 'jm';'wavelet-fft';
ring_filter_waveletfft_dec_levels = 2:6; % decomposition levels for 'wavelet-fft'
ring_filter_waveletfft_wname = 'db30';'db25'; % wavelet type for 'wavelet-fft'
ring_filter_waveletfft_sigma = 2.4; %  suppression factor for 'wavelet-fft'
ring_filter_jm_median_width = 11; % [3 11 21 31 39];
% PHASE RETRIEVAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
do_phase_retrieval = 0; % See 'PhaseFilter' for detailed description of parameters !
phase_retrieval_before = 1; % before stitching, interactive mode, etc.
phase_bin = 1; % Binning factor after phase retrieval, but before tomographic reconstruction
phase_retrieval_method = 'tie';'qpcut';  'tie'; %'qp' 'ctf' 'tie' 'qp2' 'qpcut'
phase_retrieval_reg_par = .5; % regularization parameter. larger values tend to blurrier images. smaller values tend to original data.
phase_retrieval_bin_filt = 0.15; % threshold for quasiparticle retrieval 'qp', 'qp2'
phase_retrieval_cutoff_frequ = 1 * pi; % in radian. frequency cutoff in Fourier space for 'qpcut' phase retrieval
phase_padding = 1; % padding of intensities before phase retrieval, 0: no padding
% TOMOGRAPHY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
do_tomo = 1; % run tomographic reconstruction
vol_shape = [];% shape of the volume to be reconstructed, either in absolute number of voxels or in relative number w.r.t. the default volume which is given by the detector width and height
vol_size = []; % if empty, unit voxel size is assumed
rot_angle_full = []; % in radians: empty ([]), full angle of rotation, or array of angles. if empty full rotation angles is determined automatically to pi or 2 pi
rot_angle_offset = pi; % global rotation of reconstructed volume
rot_axis_offset = []; % if empty use automatic computation
rot_axis_pos = []; % if empty use automatic computation. either offset or pos has to be empty. can't use both
rot_corr_area1 = []; % ROI to correlate projections at angles 0 & pi. Use [0.75 1] or so for scans with an excentric rotation axis
rot_corr_area2 = []; % ROI to correlate projections at angles 0 & pi
rot_corr_gradient = 0; % use gradient of intensity maps if signal variations are too weak to correlate projections
rot_axis_tilt = 0; % in rad. camera tilt w.r.t rotation axis. if empty calculate from registration of projections at 0 and pi
fbp_filter_type = 'Ram-Lak';'linear'; % Ram-Lak according to Kak/Slaney
fpb_filter_freq_cutoff = 1; % Cut-off frequency in Fourier space of the above FBP filter
fbp_filter_padding = 1; % symmetric padding for consistent boundary conditions, 0: no padding
fbp_filter_padding_method = 'symmetric';
butterworth_filter = 1; % use butterworth filter in addition to FBP filter
butterworth_order = 1;
butterworth_cutoff_frequ = 0.9;
astra_pixel_size = 1; % size of a detector pixel: if different from one 'vol_size' needs to be ajusted
take_neg_log = []; % take negative logarithm. if empty, use 1 for attenuation contrast, 0 for phase contrast
% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out_path = '';% absolute path were output data will be stored. !!overwrites the write_to_scratch flag. if empty uses the beamtime directory and either 'processed' or 'scratch_cc'
write_to_scratch = 0; % write to 'scratch_cc' instead of 'processed'
write_flatcor = 0; % save preprocessed flat corrected projections
write_phase_map = 0; % save phase maps (if phase retrieval is not 0)
write_sino = 0; % save sinograms (after preprocessing & before FBP filtering and phase retrieval)
write_sino_phase = 0; % save sinograms of phase maps
write_reco = 1; % save reconstructed slices (if do_tomo=1)
write_float = 1; % single precision (32-bit float) tiff
write_16bit = 0;
write_8bit = 0;
reco_bin = 2; % binning factor of reconstructed volume
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
parfolder = '';sprintf( 'cor_%s', correlation_method);''; % parent folder for 'reco', 'sino', 'phase', and 'flat_corrected'
subfolder_flatcor = ''; % subfolder in 'flat_corrected'
subfolder_phase_map = ''; % subfolder in 'phase_map'
subfolder_sino = ''; % subfolder in 'sino'
subfolder_reco = ''; % subfolder in 'reco'
% INTERACTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
verbose = 1; % print information to standard output
visual_output = 1; % show images and plots during reconstruction
interactive_determination_of_rot_axis = 1; % reconstruct slices with different rotation axis offsets
interactive_determination_of_rot_axis_tilt = 1; % reconstruct slices with different offset AND tilts of the rotation axis
lamino = 0; % find laminography tilt instead camera rotation
fixed_tilt = 0; % fixed other tilt
interactive_determination_of_rot_axis_slice = 0.5; % slice number, default: 0.5. if in [0,1): relative, if in (1, N]: absolute
% HARDWARE / SOFTWARE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
poolsize = 0.80; % number of workers used in a parallel pool. if > 1: absolute number. if 0 < poolsize < 1: relative amount of all cores to be used
link_data = 1; % ASTRA data objects become references to Matlab arrays. Reduces memory issues.
gpu_index = []; % GPU Device index to use, Matlab notation: index starts from 1. default: [], uses all

%% END OF PARAMETERS / SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% External call: parameters set by 'p05_reco_loop' %%%%%%%%%%%%%%%%%%%%%%%
if exist( 'external_parameter' ,'var')
    visual_output = 0;
    interactive_determination_of_rot_axis = 0;
    fields = fieldnames( external_parameter );
    for nn = 1:numel( fields )
        var_name = fields{nn};
        var_val = getfield( external_parameter, var_name );
        assignin('caller', var_name, var_val )
    end
    clear external_parameter;
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
PrintVerbose(verbose, 'Start reconstruction of ')
warning( 'off', 'MATLAB:imagesci:rtifc:missingPhotometricTag');
cdscandir = cd( scan_path );
NameCellToMat = @(name_cell) reshape(cell2mat(name_cell), [numel(name_cell{1}), numel(name_cell)])';
imsc1 = @(im) imsc( flipud( im' ) );
astra_clear % if reco was aborted, ASTRA memory is not cleared
if ~isempty( rot_axis_offset ) && ~isempty( rot_axis_pos )
    error('rot_axis_offset (%f) and rot_axis_pos (%f) cannot be used simultaneously. One must be empty.', rot_axis_offset, rot_axis_pos)
end
if abs(excentric_rot_axis)
    if crop_at_rot_axis(1) && stitch_projections(1)
        error( 'Do not use ''stitch projections'' in combination with ''crop_at_rot_axis''. Cropping at rot axis only makes sense without stitching in order to avoid oversampling artefacts within the overlap region.' )
    end
end

%% Folders: input
while scan_path(end) == filesep
    scan_path(end) = [];
end
[raw_path, scan_name] = fileparts(scan_path);
% Save raw path in file
filename = [userpath, filesep, 'experiments/p05/path_to_latest_raw'];
fid = fopen( filename , 'w' );
fprintf( fid, '%s', raw_path );
fclose( fid );
% Scan path
scan_path = [scan_path, filesep];
[beamtime_path, raw_folder] = fileparts(raw_path);
[~, beamtime_id] = fileparts(beamtime_path);
if ~strcmp(raw_folder, 'raw')
    error('Name of folder is not raw: %s', raw_folder)
end
PrintVerbose(verbose, '%s', scan_name)
PrintVerbose(verbose, '\n scan_path:%s', scan_path)
% Save scan path to file
filename = [userpath, filesep, 'experiments/p05/path_to_latest_scan'];
fid = fopen( filename , 'w' );
fprintf( fid, '%s', scan_path );
fclose( fid );

%% Folders: output
out_folder = 'processed';
if write_to_scratch(1)
    out_folder = 'scratch_cc';
end
if isempty( out_path )
    out_path = [beamtime_path, filesep, out_folder, filesep, scan_name];
else
    out_path = [out_path, filesep, scan_name];
end
if ~isempty(parfolder)
    out_path = [out_path, filesep, parfolder];
end

% Projections path
if isempty( subfolder_flatcor )
    flatcor_path = [out_path, filesep, 'flat_corrected', filesep];
else
    flatcor_path = [out_path, filesep, 'flat_corrected', filesep, subfolder_flatcor, filesep];
end
PrintVerbose(verbose & write_flatcor, '\n flatcor_path: %s', flatcor_path)

% Phase map path
if isempty( subfolder_phase_map )
    phase_map_path = [out_path, filesep, 'phase_map', filesep];
else
    phase_map_path = [out_path, filesep, 'phase_map', filesep, subfolder_phase_map, filesep];
end
PrintVerbose(verbose & write_phase_map, '\n phase_map_path: %s', phase_map_path)

% Sinogram path
if isempty( subfolder_sino )
    sino_path = [out_path, filesep, 'sino', filesep];
    sino_phase_path = [out_path, filesep, 'sino_phase', filesep];
else
    sino_path = [out_path, filesep, 'sino', filesep, subfolder_sino, filesep];
    sino_phase_path = [out_path, filesep, 'sino_phase', filesep, subfolder_sino, filesep];
end
PrintVerbose(verbose & write_sino, '\n sino_path: %s', sino_path)
PrintVerbose(verbose & write_sino_phase, '\n sino_phase_path: %s', sino_phase_path)

% Reco path
if isempty( subfolder_reco )
    reco_path = [out_path, filesep, 'reco', filesep];
else
    reco_path = [out_path, filesep, 'reco', filesep, subfolder_reco, filesep];
end

% Memory
PrintVerbose(verbose, '\n hostname : %s', getenv( 'HOSTNAME' ) );
[mem_free, mem_avail, mem_total] = free_memory;
PrintVerbose(verbose, '\n system memory: free, available, total : %.3g GiB, %.3g GiB, %.3g GiB', mem_free/1024^3, mem_avail/1024^3, mem_total/1024^3)
if isempty( gpu_index )
    gpu_index = 1:gpuDeviceCount;
end
for nn = gpu_index
    gpu = parallel.gpu.GPUDevice.getDevice( nn );
    PrintVerbose(verbose, '\n gpu %u memory : available, total, percent = %.3g GiB, %.3g GiB, %.2f%%', nn, gpu.AvailableMemory/1024^3, gpu.TotalMemory/1024^3, 100*gpu.AvailableMemory/gpu.TotalMemory)
end

%% File names

% Projection file names
proj_names = FilenameCell( [scan_path, '*.img'] );
if isempty( proj_names )
    proj_names =  FilenameCell( [scan_path, 'proj_*.tif'] );
end
if isempty( proj_names )
    proj_names =  FilenameCell( [scan_path, '*img*.raw'] );
end
num_proj_found = numel(proj_names);

% Ref file names
ref_names = FilenameCell( [scan_path, '*.ref'] );
if isempty( ref_names )
    ref_names = FilenameCell( [scan_path, 'ref_*.tif'] );
end
if isempty( ref_names )
    ref_names =  FilenameCell( [scan_path, '*ref*.raw'] );
end
num_ref_found = numel(ref_names);
if isempty( ref_range )
    ref_range = 1;
end
if numel( ref_range ) == 1
    ref_range = 1:ref_range:num_ref_found;
end
num_ref_used = numel( ref_range );
ref_names_mat = NameCellToMat( ref_names(ref_range) );
ref_nums = CellString2Vec( ref_names(ref_range) );
PrintVerbose(verbose, '\n number of refs found : %g', num_ref_found)
PrintVerbose(verbose, '\n number of refs used : %g', num_ref_used)
PrintVerbose(verbose, '\n reference range used : %g:%g:%g%', ref_range(1), ref_range(2) - ref_range(1), ref_range(end))

% Dark file names
dark_names = FilenameCell( [scan_path, '*.dar'] );
if isempty( dark_names )
    dark_names = FilenameCell( [scan_path, 'dark_*.tif'] );
end
if isempty( dark_names )
    dark_names =  FilenameCell( [scan_path, '*dar*.raw'] );
end
dark_nums = CellString2Vec( dark_names );
num_dark = numel(dark_names);
PrintVerbose(verbose, '\n number of darks found : %g', num_dark)

% Projection range to read
if isempty( proj_range )
    proj_range = 1;
end
if numel( proj_range ) == 1
    proj_range = 1:proj_range:num_proj_found;
end
num_proj_used = numel( proj_range );
proj_nums = CellString2Vec( proj_names(proj_range) );
PrintVerbose(verbose, '\n number of projections found : %g', num_proj_found)
PrintVerbose(verbose, '\n number of projections used : %g', num_proj_used)
PrintVerbose(verbose, '\n projection range used : first:stride:last =  %g:%g:%g', proj_range(1), proj_range(2) - proj_range(1), proj_range(end))

h5log = sprintf('%s%s_nexus.h5',scan_path,scan_name);
raw_im_shape_uncropped = [];
dtype = '';
tif_info = [];
offset_shift = 0;
if ~exist( h5log, 'file')
    
    %% Image shape and ROI
    filename = sprintf('%s%s', scan_path, ref_names{1});
    [im_raw, tif_info] = read_image( filename );
    raw_im_shape_uncropped = size( im_raw );
    im_roi = read_image( filename, '', raw_roi, tif_info );
    raw_im_shape = size( im_roi );
    raw_im_shape_binned = floor( raw_im_shape / raw_bin );
    raw_im_shape_binned1 = raw_im_shape_binned(1);
    raw_im_shape_binned2 = raw_im_shape_binned(2);
    
    %% P05 log-file
    str = dir( sprintf( '%s*scan.log', scan_path) );
    filename = sprintf( '%s/%s', str.folder, str.name);
    [par, cur, cam] = p05_log( filename );
    if isempty( energy )
        energy = par.energy;
    end
    if isempty( eff_pixel_size )
        eff_pixel_size = par.eff_pixel_size;
    end
    eff_pixel_size_binned = raw_bin * eff_pixel_size;
    if isempty( sample_detector_distance )
        sample_detector_distance = par.sample_detector_distance;
    end
else
    
    %% P05 standard log file
    str = dir( sprintf( '%s*scan.log', scan_path) );
    filename = sprintf( '%s/%s', str.folder, str.name);
    [par, cur, cam] = p05_log( filename );
    eff_pixel_size = par.eff_pixel_size;
    eff_pixel_size_binned = raw_bin * eff_pixel_size;
    sample_detector_distance = par.o_ccd_dist / 1000;
    
    %% HDF5 log
    h5i = h5info( h5log );
    
    clear stimg_name stimg_key petra s_rot s_stage_x
    
    stimg_name.value = h5read( h5log, '/entry/scan/data/image_file/value');
    stimg_name.time = h5read( h5log,'/entry/scan/data/image_file/time');
    stimg_key.value = h5read( h5log,'/entry/scan/data/image_key/value');
    stimg_key.time = h5read( h5log,'/entry/scan/data/image_key/time');
    
    petra.time = h5read( h5log,'/entry/hardware/beam_current/current/time');
    petra.current = h5read( h5log,'/entry/hardware/beam_current/current/value');
    
    s_rot.time = h5read( h5log, '/entry/scan/data/s_rot/time');
    s_rot.value = h5read( h5log, '/entry/scan/data/s_rot/value');
    
    s_stage_x.time = h5read( h5log, '/entry/scan/data/s_stage_x/time');
    s_stage_x.value = h5read( h5log, '/entry/scan/data/s_stage_x/value');
    
    offset_shift = s_stage_x.value( ~boolean( stimg_key.value(par.n_dark+1:end) ) ) * 1e-3 / eff_pixel_size_binned;
    
    %% Image shape
    if strcmp( cam, 'EHD')
        raw_im_shape_uncropped = [3056 3056];
        dtype = 'uint16';
    end
    raw_im_shape = raw_im_shape_uncropped;
    filename = sprintf('%s%s', scan_path, ref_names{1});
    im_raw = read_raw( filename, raw_im_shape, 'uint16' );
    raw_im_shape_binned = floor( raw_im_shape / raw_bin );
    raw_im_shape_binned1 = raw_im_shape_binned(1);
    raw_im_shape_binned2 = raw_im_shape_binned(2);
    
    ring_filter = 0; % % not needed
    
    %% TODO: FIX
    energy = 1; % NOT LOGGED AT ALL
    ring_current_normalization = 0;
    rot_angle_full = pi; % FIX
    
end
PrintVerbose( verbose, '\n energy : %.1f keV', energy / 1e3 )
PrintVerbose( verbose, '\n distance sample dector : %.1f mm', sample_detector_distance * 1000 )
PrintVerbose( verbose, '\n effective pixel size unbinned : %.2f micron',  eff_pixel_size * 1e6)
PrintVerbose(verbose, '\n raw image shape : %g  %g', raw_im_shape_uncropped)
PrintVerbose(verbose, '\n raw image shape roi : %g  %g', raw_im_shape)
PrintVerbose(verbose, '\n raw image shape binned : %g  %g', raw_im_shape_binned)
PrintVerbose(verbose, '\n raw binning factor : %u', raw_bin)
PrintVerbose(verbose, '\n raw bining before pixel filtering : %u', bin_before_filtering)

%% Start parallel CPU pool
t = toc;
if poolsize <= 1 && poolsize > 0
    poolsize = max( floor( poolsize * feature('numCores') ), 1 );
end
PrintVerbose( verbose && (poolsize > 1), '\nStart parallel pool of %u workers. ', poolsize)
OpenParpool(poolsize);
PrintVerbose( verbose && (poolsize > 1), ' Time elapsed: %.1f s', toc-t)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read flat corrected projection
if read_flatcor(1)
    t = toc;
    if isempty( read_flatcor_path )
        read_flatcor_path = flatcor_path;
    end
    % File names
    data_struct = dir( [read_flatcor_path filesep 'proj*.*'] );
    if isempty( data_struct )
        fprintf('\n No flat corrected projections found! Switch to standard pre-processing.')
        read_flatcor = 0;
    else
        proj_names = {data_struct.name};
        num_proj_read = numel(proj_names);
        if num_proj_used > num_proj_read
            fprintf('\n Less projections available (%g) than demanded (%g)! Switch to standard pre-processing.', num_proj_read, num_proj_used )
            read_flatcor = 0;
        end
    end
    PrintVerbose(verbose, '\n Read flat corrected projections.')
    if num_proj_read ~= num_proj_used
        fprintf('\n WARNING: Number of flat corrected projections read (%g) differs from number of projections to be processed (%g)!\n', num_proj_read, num_proj_used)
    end
    
    % Preallocation
    proj = zeros( raw_im_shape_binned(1), raw_im_shape_binned(2), num_proj_read, 'single');
    PrintVerbose( verbose, ' Allocated bytes: %.2f GiB.', Bytes( proj, 3 ) )
    proj_names_mat = NameCellToMat( proj_names );
    
    % Read flat corrected projections
    parfor nn = 1:num_proj_read
        filename = sprintf('%s%s', flatcor_path, proj_names_mat(nn, :));
        proj(:, :, nn) = fliplr( read_image( filename ) );
    end
    PrintVerbose(verbose, ' Time elapsed: %.1f s (%.2f min)', toc - t, ( toc - t ) / 60 )
    
    %% Read raw data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ~read_flatcor
    %% Dark field
    t = toc;
    PrintVerbose(verbose, '\nProcessing dark fields.')
    darks = zeros( [raw_im_shape_binned, num_dark], 'single');
    PrintVerbose( verbose, ' Allocated bytes: %.2f MiB.', Bytes( darks, 2 ) )
    parfor nn = 1:num_dark
        filename = sprintf('%s%s', scan_path, dark_names{nn});
        im = single( read_image( filename, '', raw_roi, tif_info, raw_im_shape_uncropped, dtype) );
        % Remove large outliers. Assume Poisson distribtion at large lambda
        % is approximately a Gaussian distribution and set all value above
        % mean + 4 * std (99.994 of values lie within 4 std). Due to
        % outliers 4*std will contain much more values and is a good
        % estimate
        im_mean = mean( im(:) );
        im_std = std( im(:) );
        im( im > im_mean + 4*im_std) = im_mean;
        if bin_before_filtering
            darks(:, :, nn) = FilterPixel( Binning(im, raw_bin) / raw_bin^2, dark_FiltPixThresh );
        else
            darks(:, :, nn) = Binning( FilterPixel( im, dark_FiltPixThresh), raw_bin) / raw_bin^2;
        end
    end
    darks_min = min( darks(:) );
    darks_max = max( darks(:) );
    % Reject dark images which are all zero
    darks_to_use = zeros( 1, num_dark, 'logical' );
    parfor nn = 1:num_dark
        darks_to_use(nn) = boolean( max2( darks(:,:,nn) )  );
    end
    % Median dark
    dark = squeeze( median(darks(:,:,darks_to_use), 3) );
    clear darks
    dark_med_min = min( dark(:) );
    dark_med_max = max( dark(:) );
    PrintVerbose(verbose, ' Time elapsed: %.1f s', toc-t)
    if visual_output(1)
        h1 = figure('Name', 'data and flat-and-dark-field correction');
        subplot(2,3,1)
        imsc1( dark );
        axis equal tight
        xticks([]);yticks([])
        title(sprintf('median dark field'))
        colorbar
        drawnow
    end
    
    %% Flat field
    t = toc;
    PrintVerbose(verbose, '\nProcessing flat fields.')
    
    % Preallocation
    flat = zeros( [raw_im_shape_binned, num_ref_used], 'single');
    PrintVerbose( verbose, ' Allocated bytes: %.2f GiB.', Bytes( flat, 3 ) )
    
    % Parallel loop
    refs_to_use = zeros( 1, size( flat,3), 'logical');
    parfor nn = 1:num_ref_used
        filename = sprintf('%s%s', scan_path, ref_names_mat(nn,:));
        if bin_before_filtering
            flat(:, :, nn) = FilterPixel( Binning( single( read_image( filename, '', raw_roi, tif_info, raw_im_shape_uncropped, dtype ) ), raw_bin) / raw_bin^2, ref_FiltPixThresh);
        else
            flat(:, :, nn) = Binning( FilterPixel( read_image( filename, '', raw_roi, tif_info, raw_im_shape_uncropped, dtype ), ref_FiltPixThresh), raw_bin) / raw_bin^2;
        end
        
        % Check for zeros
        num_zeros =  sum( sum( flat(:,:,nn) < 1 ) );
        refs_to_use(nn) = ~boolean( num_zeros  );
        %         if num_zeros > 0
        %             fprintf('\n WARNING: values of %u pixels of flat field no %u are below 1.\n', num_zeros, nn)
        %         end
    end
    
    % Delete empty refs
    flat(:,:,~refs_to_use) = [];
    
    % min/max values before dark field subtraction and ring current normalization
    flat_min = min( flat(:) );
    flat_max = max( flat(:) );
    
    % Dark field correction
    flat = bsxfun( @minus, flat, dark );
    
    % Ring current normalization
    if ring_current_normalization(1)
        switch lower( cam )
            case 'kit'
                %ref_ind = ref_nums + 1 ;
                ref_check = ref_range - 1;
            case 'ehd'
                ref_check = ref_nums;
        end
        if isequal( ref_check, [cur.ref(ref_range).ind])
            ref_rc = [cur.ref(ref_range).val];
            ref_rcm = mean( ref_rc(:) );
            scale_factor = 100 ./ shiftdim( ref_rc(refs_to_use), -1 );
            flat = bsxfun( @times, flat, scale_factor );
            if visual_output(1)
                hrc = figure('Name', 'Ring currents');
                subplot(2,1,1);
                plot( ref_rc(:), '.' )
                axis tight
                title(sprintf('ring current: flat fields'))
                legend( sprintf( 'mean: %.2f mA', ref_rcm) )
                drawnow
            end
        else
            fprintf('\n WARNING: flat fields not normalized by ring current. Names read from dir() and log-file are inconsistent.\n')
        end
    end
    
    % Enforce positivity
    parfor nn = 1:size( flat, 3 )
        im = flat(:,:,nn);
        m = im < 1;
        im(m) = 1;
        flat(:,:,nn) = im;
    end
    flat_min2 = min( flat(:) );
    flat_max2 = max( flat(:) );
    num_zeros =  sum( flat(:) < 1 );
    if num_zeros > 0
        fprintf('\n WARNING: flat field contains %u zeros\n', num_zeros)
    end
    
    nn = sum( ~refs_to_use(:) );
    num_ref_used = num_ref_used - nn;
    PrintVerbose(verbose, ' Time elapsed: %.1f s', toc-t)
    PrintVerbose( verbose && nn,'\n discarded empty refs : %u', nn )
    
    % Show flat field
    if visual_output(1)
        if exist( 'h1' , 'var' )
            figure(h1)
        else
            h1 = figure('Name', 'data and flat-and-dark-field correction');
        end
        subplot(2,3,2)
        imsc1( flat(:,:,1) )
        axis equal tight
        xticks([]);yticks([])
        title(sprintf('flat field #1'))
        drawnow
        colorbar
    end
    
    %% Projections
    t = toc;
    PrintVerbose(verbose, '\nRead and filter projections.')
    % Preallocation
    proj = zeros( raw_im_shape_binned(1), raw_im_shape_binned(2), num_proj_used, 'single');
    PrintVerbose( verbose, ' Allocated bytes: %.2f GiB.', Bytes( proj, 3 ) )
    img_names_mat = NameCellToMat( proj_names(proj_range) );
    
    % Display first raw image
    if visual_output(1)
        if exist( 'h1' , 'var' )
            figure(h1)
        else
            h1 = figure('Name', 'data and flat-and-dark-field correction');
        end
        filename = sprintf('%s%s', scan_path, img_names_mat(1, :));
        raw1 = Binning( FilterPixel( read_image( filename, '', raw_roi, tif_info, raw_im_shape_uncropped, dtype ), proj_FiltPixThresh), raw_bin) / raw_bin^2;
        subplot(2,3,3)
        imsc1( raw1 )
        axis equal tight
        xticks([]);yticks([])
        title(sprintf('raw proj #1'))
        drawnow
        colorbar
    end
    
    % Get absolute filter thresholds from percentage-wise pixel filtering
    % of 1st, middle, and last projection to speed up processing
    if proj_FiltPixThresh(1) < 1 || proj_FiltPixThresh(2) < 0.5
        
        filename = sprintf('%s%s', scan_path, img_names_mat(num_proj_used, :));
        [~, ht(3), dt(3)] = FilterPixel( read_image( filename, '', raw_roi, tif_info, raw_im_shape_uncropped, dtype ), proj_FiltPixThresh);
        
        filename = sprintf('%s%s', scan_path, img_names_mat(1, :));
        [~, ht(2), dt(2)] = FilterPixel( read_image( filename, '', raw_roi, tif_info, raw_im_shape_uncropped, dtype ), proj_FiltPixThresh);
        
        filename = sprintf('%s%s', scan_path, img_names_mat(round(num_proj_used/2), :));
        [~, ht(1), dt(1)] = FilterPixel( read_image( filename, '', raw_roi, tif_info, raw_im_shape_uncropped, dtype ), proj_FiltPixThresh);
        
        HotThresh = median( ht );
        DarkThresh = median( dt );
    else
        HotThresh = proj_FiltPixThresh(1);
        DarkThresh = proj_FiltPixThresh(2);
    end
    
    % Read raw projections
    projs_to_use = zeros( 1, size( proj,3), 'logical' );
    parfor nn = 1:num_proj_used
        % Read and filter raw projection
        filename = sprintf('%s%s', scan_path, img_names_mat(nn,:));
        if bin_before_filtering
            im = FilterPixel( Binning( single( read_image( filename, '', raw_roi, tif_info, raw_im_shape_uncropped, dtype ) ), raw_bin) / raw_bin^2, [HotThresh, DarkThresh]);
        else
            im = Binning( FilterPixel( read_image( filename, '', raw_roi, tif_info, raw_im_shape_uncropped, dtype ), [HotThresh, DarkThresh]), raw_bin) / raw_bin^2;
        end
        
        % Check for zeros and reject images which is all zero
        num_zeros =  sum( sum( im < 1 ) );
        projs_to_use(nn) = ~boolean( num_zeros  );
        
        if projs_to_use(nn)
            proj(:, :, nn) = im;
        end
    end
    
    % Delete empty projections
    proj(:,:,~projs_to_use) = [];
    
    raw_min = min( proj(:) );
    raw_max = max( proj(:) );
    
    % Dark field correction
    proj = bsxfun( @minus, proj, dark);
    
    % Ring current normalization
    if ring_current_normalization(1)
        switch lower( cam )
            case 'kit'
                %proj_ind = proj_nums + 1;
                proj_check = proj_range - 1;
            case 'ehd'
                proj_check = proj_nums;
        end
        if isequal( proj_check, [cur.proj(proj_range).ind] )
            proj_rc = [cur.proj(proj_range).val];
            proj_rcm = mean( proj_rc(:) );
            scale_factor = 100 ./ shiftdim( proj_rc(projs_to_use), -1 );
            proj = bsxfun( @times, proj, scale_factor );
            % Plot ring current
            if visual_output(1)
                if exist( 'hrc', 'var' )
                    figure(hrc)
                else
                    hrc = figure('Name', 'Ring currents');
                end
                subplot(2,1,2);
                plot( proj_rc(:), '.' )
                axis tight
                title(sprintf('ring current: raw projections'))
                legend( sprintf( 'mean: %.2f mA', proj_rcm) )
                drawnow
            end
        else
            fprintf('\n WARNING: projections not normalized by ring current. Names read from dir() and log-file are inconsistent.\n')
        end
    end
    
    % Enforce positivity
    parfor nn = 1:size( proj, 3 )
        im = proj(:,:,nn);
        m = im < 1;
        if sum( m(:) ) > 0
            im(m) = 1;
            proj(:,:,nn) = im;
        end
    end
    
    raw_min2 = min( proj(:) );
    raw_max2 = max( proj(:) );
    nn = sum( ~projs_to_use(:) );
    num_proj_used = num_proj_used - nn;
    PrintVerbose( verbose, ' Time elapsed: %.1f s (%.2f min)', toc - t, ( toc - t ) / 60 )
    PrintVerbose( verbose, '\n hot- / dark-pixel filter threshold : %f, %f', HotThresh, DarkThresh )
    PrintVerbose( verbose && nn,'\n discarded empty projections : %u', nn )
    
    %% Flat field correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STOP HERE FOR FLATFIELD CORRELATION MAPPING
    [proj, corr, proj_roi, flat_roi] = proj_flat_correlation( proj, flat, correlation_method, flat_corr_area1, flat_corr_area2, raw_im_shape_binned, corr_shift_max_pixelshift, corr_num_flats, decimal_round_precision, flatcor_path, verbose, visual_output);
    
    PrintVerbose(verbose, '\n sinogram size = [%g, %g, %g]', size( proj ) )
    if visual_output(1)
        if exist( 'h1' , 'var' )
            figure(h1)
        else
            h1 = figure('Name', 'data and flat-and-dark-field correction');
        end
        
        subplot(2,3,4)
        imsc1( proj(:,:,1))
        axis equal tight
        xticks([]);yticks([])
        title(sprintf('intensity: first proj'))
        colorbar
        
        subplot(2,3,5)
        imsc1( proj(:,:,round(size(proj,3)/2)))
        axis equal tight
        xticks([]);yticks([])
        title(sprintf('intensity: middle proj'))
        colorbar
        
        subplot(2,3,6)
        imsc1( proj(:,:,end))
        axis equal tight
        xticks([]);yticks([])
        title(sprintf('intensity: last proj'))
        colorbar
        
        drawnow
    end
    
    %% Angles
    if exist('cur', 'var') && isfield(cur, 'proj') && isfield( cur.proj, 'angle')
        % for KIT cam this includes missing angles
        angles = [cur.proj.angle] / 180 * pi;
        if strcmpi(cam, 'kit')
            % drop angles where projections are missing
            angles = angles(1 + proj_nums);
        else
            angles = angles(proj_range);
        end
    elseif exist( h5log, 'file')
        num_proj = double( par.n_angle );
        angles = s_rot.value( ~boolean( stimg_key.value(par.n_dark+1:end) ) ) * pi / 180;
    else
        if isfield( par, 'num_projections' )
            num_proj = double( par.num_projections );
        elseif isfield( par, 'projections' )
            num_proj = double( par.projections );
        elseif isfield( par, 'n_angles' )
            num_proj = double( par.n_angles );
        end
        switch lower( cam )
            case 'ehd'
                angles = rot_angle_full * (0:num_proj - 1) / (num_proj - 1); % EHD: ok
            case 'kit'
                angles = rot_angle_full * (0:num_proj - 1) / num_proj; % KIT: ok if par.projections exist
        end
    end
    % drop angles where projections are empty
    angles(~projs_to_use) = [];
    % sort angles for displaying and non-contiguous acquisition
    [~, sorted_angle_index] = sort( angles );
    
    %% Ring artifact filter
    sino_slice = round( size( proj, 2) / 2 );
    sino_unfilt = squeeze( proj(:,sino_slice,sorted_angle_index) )';
    if ring_filter(1)
        t = toc;
        PrintVerbose(verbose, '\nFilter ring artifacts using %s. ', ring_filter_method)
        switch lower( ring_filter_method )
            
            % Combined wavelet-FFT filter to remove ring artrifacts
            case 'wavelet-fft'
                
                parfor nn = 1:size( proj, 2)
                    sino = squeeze( proj(:,nn,:) )';
                    [d2, d1] = size( sino );
                    sino = FilterStripesCombinedWaveletFFT( sino, ring_filter_waveletfft_dec_levels, ring_filter_waveletfft_wname, ring_filter_waveletfft_sigma )';
                    sino = sino( 1:d1, 1:d2 );
                    proj(:,nn,:) = sino;
                end
                PrintVerbose(verbose, ' Time elapsed: %.1f s (%.2f min)', toc-t, (toc-t)/60)
                sino_filt = squeeze( proj(:,sino_slice,sorted_angle_index) )';
                
                if visual_output(1)
                    
                    h3 = figure('Name', 'Sinogram and ring filter');
                    
                    subplot(3,1,1)
                    imsc( sino_unfilt )
                    axis equal tight
                    title(sprintf('sino unfiltered, y = %u', sino_slice))
                    colorbar
                    
                    subplot(3,1,2)
                    imsc( sino_filt )
                    axis equal tight
                    title(sprintf('sino filtered, y = %u', sino_slice))
                    colorbar
                    
                    subplot(3,1,3)
                    imsc( sino_filt - sino_unfilt )
                    axis equal tight
                    title(sprintf('sino filt - sino, y = %u', sino_slice))
                    colorbar
                    
                    drawnow
                end
                
                % Simple ring artifact filter
            case 'jm'
                if numel( ring_filter_jm_median_width ) > 1
                    %% Combine if/else
                    for nn = ring_filter_jm_median_width
                        %m=sum(im);
                        proj_mean = sum( proj, 3);
                        
                        %md=medfilt1(m,nn);
                        proj_mean_med = medfilt2( proj_mean, [nn, 1], 'symmetric' );
                        
                        %f=md./m;
                        mask = proj_mean_med ./ proj_mean;
                        
                        %im=bsxfun(@times,im,f);
                        proj = bsxfun( @times, proj, mask);
                    end
                else
                    proj_mean = mean( proj, 3);
                    proj_mean_med = medfilt2( proj_mean, [ring_filter_jm_median_width, 1], 'symmetric' );
                    mask = proj_mean_med ./ proj_mean;
                    proj = bsxfun( @times, proj, mask);
                end
                sino_filt = squeeze( proj(:,sino_slice,sorted_angle_index) )';
                
                PrintVerbose(verbose, ' Time elapsed: %.1f s (%.2f min)', toc-t, (toc-t)/60)
                PrintVerbose( verbose, '\n ring filter mask min/max: %f, %f', min( mask(:) ), max( mask(:) ) )
                if visual_output(1)
                    h3 = figure('Name', 'Sinogram and ring filter');
                    
                    subplot(2,2,1)
                    imsc( sino_unfilt )
                    axis equal tight
                    title(sprintf('sino unfiltered, y = %u', sino_slice))
                    colorbar %('Location', 'southoutside')
                    
                    subplot(2,2,2)
                    imsc( sino_filt - sino_unfilt )
                    axis equal tight
                    title(sprintf('sino filt - sino, y = %u', sino_slice))
                    colorbar %('Location', 'southoutside')
                    
                    subplot(2,2,3)
                    imsc1( proj_mean )
                    axis equal tight
                    title(sprintf('mean projection'))
                    colorbar %('Location', 'southoutside')
                    
                    subplot(2,2,4)
                    imsc1( FilterHisto( mask, 5 ) )
                    axis equal tight
                    title(sprintf('mask for normalization'))
                    colorbar %('Location', 'southoutside')
                    
                    drawnow
                end
        end
    end
    proj_min = min( proj(:) );
    proj_max = max( proj(:) );
    
    %% Write corrected projections
    if write_flatcor(1)
        t = toc;
        PrintVerbose(verbose, '\nSave flat-corrected projections.')
        CheckAndMakePath( flatcor_path )
        % Projections
        parfor nn = 1:size( proj, 3 )
            filename = sprintf('%sproj_%06u.tif', flatcor_path, nn );
            write32bitTIFfromSingle(filename, rot90( proj(:, :, nn) ) );
        end
        PrintVerbose(verbose, ' Time elapsed: %.1f (%.2f min)', toc-t, (toc-t)/60)
    end
end

%% Rotation axis position and tomgraphic reconstruction parameters %%%%%%%%
if do_phase_retrieval(1)
    if isempty( take_neg_log )
        take_neg_log = 0;
    end
else
    phase_bin = 0;
end

%% Phase retrieval %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
edp = [energy, sample_detector_distance, eff_pixel_size_binned];
if do_phase_retrieval(1) && phase_retrieval_before(1)
    PrintVerbose(verbose, '\nPhase retrieval')
    t = toc;
    PrintVerbose(verbose, '\n energy : %g eV', energy)
    PrintVerbose(verbose, '\n sample detector distance : %g m', sample_detector_distance)
    
    if isempty( take_neg_log )
        take_neg_log = 0;
    end
    
    % Phase retrieval filter
    im_shape = [size(proj,1) , size(proj,2)];
    im_shape_pad = (1 + phase_padding) * im_shape;
    [phase_filter, pha_appendix] = PhaseFilter( phase_retrieval_method, im_shape_pad, edp, phase_retrieval_reg_par, phase_retrieval_bin_filt, phase_retrieval_cutoff_frequ, 'single');
    
    % reco phase dir
    if isempty( subfolder_reco )
        reco_phase_path = [out_path, filesep, 'reco_phase', filesep, pha_appendix, filesep];
    else
        reco_phase_path = [out_path, filesep, 'reco_phase', filesep, pha_appendix, filesep, subfolder_reco, filesep];
    end
    CheckAndMakePath( reco_phase_path )
    PrintVerbose(verbose, '\n reco_phase_path : %s', reco_phase_path)
    PrintVerbose(verbose, '\n phase retrieval method : %s', phase_retrieval_method)
    
    % Retrieval
    parfor nn = 1:size(proj, 3)
        % combined GPU and parfor usage requires memory management
        im = padarray( proj(:,:,nn), phase_padding * im_shape, 'symmetric', 'post' );
        %im = padarray( gpuArray( proj(:,:,nn) ), raw_im_shape_binned, 'post', 'symmetric' );
        im = -real( ifft2( phase_filter .* fft2( im ) ) );
        proj(:,:,nn) = im(1:im_shape(1), 1:im_shape(2));
        %proj(:,:,nn) = gather( im(1:raw_im_shape_binned1, 1:raw_im_shape_binned2) );
    end
    pause(0.01)
    PrintVerbose(verbose, '\n done in %g s (%.2f min)', toc-t, (toc-t)/60)
    
    % Save phase maps
    if write_phase_map(1)
        t = toc;
        PrintVerbose(verbose, '\nSave phase maps:')
        phase_map_path = [phase_map_path, pha_appendix, filesep];
        CheckAndMakePath( phase_map_path );
        parfor nn = 1:size( proj, 3)
            filename = sprintf( '%sphase_%06u.tif', phase_map_path, nn);
            write32bitTIFfromSingle( filename, squeeze( rot90( proj( :, :, nn) ) ) )
        end
        pause(0.01)
        PrintVerbose(verbose, ' Time elapsed: %.1f s (%.2f min)', toc-t, (toc-t)/60)
    end
end

tint = 0;
if do_tomo(1)
    t = toc;
    
    % ROI for correlation of projections at angles 0 & pi
    if isempty( rot_corr_area1 )
        switch excentric_rot_axis
            case -1
                rot_corr_area1 = [0 0.25];
            case 0
                rot_corr_area1 = [0.25 0.75];
            case 1
                rot_corr_area1 = [0.75 1];
        end
    end
    if isempty( rot_corr_area2 )
        rot_corr_area2 = [0.1 0.9];
    end
    rot_corr_area1 = IndexParameterToRange( rot_corr_area1, raw_im_shape_binned1 );
    rot_corr_area2 = IndexParameterToRange( rot_corr_area2, raw_im_shape_binned2 );
    
    % Full rotation angle
    if isempty( rot_angle_full )
        if isfield( par, 'rotation')
            % From log file
            rot_angle_full = par.rotation / 180 * pi;
        elseif exist('cur', 'var') && isfield(cur, 'proj') && isfield( cur.proj, 'angle')
            % from beam current log
            rot_angle_full = (cur.proj(end).angle - cur.proj(1).angle) * pi /180; % KIT: , EHD: ok
        else
            % Guess from correlation of projections
            im1 = proj( rot_corr_area1, rot_corr_area2, 1);
            im2 = proj( rot_corr_area1, rot_corr_area2, end);
            [~, cm1] = ImageCorrelation( im1, im2, 0, 0, 0); % large if full angle is  2 * pi
            [~, cm2] = ImageCorrelation( im1, flipud(im2), 0, 0, 0); % large if full angle is pi
            if max( cm1(:) ) > max( cm2(:) )
                rot_angle_full = 2 * pi;
            elseif max( cm1(:) ) < max( cm2(:) )
                rot_angle_full = pi;
            else
                error('Determination of full angle of rotation from correlation of first and last projection not successful.')
            end
        end
    end
    PrintVerbose(verbose, '\n full rotation angle: %g * pi', rot_angle_full / pi)
    if numel( angles ) ~= num_proj_used
        error('Number of elements in array of angles (%g) unequal number of projections read (%g)', numel( angles ), num_proj_used)
    end
    % retrieve index at angles 0 and pi
    [val1, ind1] = min( abs( angles ));
    [val2, ind2] = min( abs( angles - pi ));
    % ROI
    im1 = proj( rot_corr_area1, rot_corr_area2, ind1);
    im2 = flipud( proj( rot_corr_area1, rot_corr_area2, ind2) );
    PrintVerbose(verbose, '\n correlation of images : [filename index, projection index, angle / pi] = [%g, %g, %g] and [%g, %g, %g]', ...
        proj_nums(ind1), ind1, val1 / pi, proj_nums(ind2), ind2, (val2 + pi) / pi )
    if rot_corr_gradient(1)
        l = 2;
        [g11, g12] = gradient(im1);
        im1 = abs(g11).^l + abs(g12).^l;
        [g21, g22] = gradient(im2);
        im2 = abs(g21).^l + abs(g22).^l;
    end
    % Correlation
    out = ImageCorrelation( im1, im2, 0, 0, 0);
    % relative shift
    rot_corr_shift_x = round( out.shift1, decimal_round_precision) + rot_corr_area1(1) + (rot_corr_area1(end) - raw_im_shape_binned1) - 1;
    rot_corr_shift_y = round( out.shift2, decimal_round_precision) + rot_corr_area2(1) + (rot_corr_area1(end) - raw_im_shape_binned1) - 1;
    PrintVerbose(verbose, '\n relative shift: %g, %g', rot_corr_shift_x, rot_corr_shift_y)
    rot_axis_offset_calc = rot_corr_shift_x / 2;
    rot_axis_pos_calc = raw_im_shape_binned1 / 2 + rot_axis_offset_calc;
    PrintVerbose(verbose, '\n calculated rotation axis offset: %.2f', rot_axis_offset_calc)
    % use calculated offset if both offset and position are empty
    if isempty( rot_axis_offset ) && isempty( rot_axis_pos )
        rot_axis_offset = rot_axis_offset_calc;
    end
    if isempty( rot_axis_pos )
        rot_axis_pos = raw_im_shape_binned1 / 2 + rot_axis_offset;
    end
    if isempty( rot_axis_offset )
        rot_axis_offset = rot_axis_pos - raw_im_shape_binned1 / 2;
    end
    
    % Tilt of rotation axis
    if interactive_determination_of_rot_axis_tilt(1)
        im1c = RotAxisSymmetricCropping( proj(:,rot_corr_area2,ind1), rot_axis_pos, 1);
        im2c = flipud(RotAxisSymmetricCropping( proj(:,rot_corr_area2,ind2) , rot_axis_pos, 1));
        [optimizer, metric] = imregconfig('monomodal');
        tform = imregtform(im2c, im1c, 'rigid', optimizer, metric);
        rot_axis_tilt_calc = asin( tform.T(1,2) ) / 2;
        PrintVerbose(verbose, '\n calculated tilt of rotation axis : %g rad = %g deg', rot_axis_tilt_calc, rot_axis_tilt_calc * 180 / pi)
    else
        rot_axis_tilt_calc = [];
    end
    if isempty( rot_axis_tilt )
        if abs( rot_axis_pos_calc ) < 0.01 % 0.573 degree
            rot_axis_tilt = rot_axis_tilt_calc;
            PrintVerbose(verbose, '\n use calculated tilt of rotation axis')
        else
            rot_axis_tilt = 0;
        end
    end
    
    %% Determine rotation axis position
    tint = 0;
    if interactive_determination_of_rot_axis(1)
        tint = toc;
        fprintf( '\n\nENTER INTERACTIVE MODE' )
        fprintf( '\n number of pixels: %u', raw_im_shape_binned1)
        fprintf( '\n image center: %.1f', raw_im_shape_binned1 / 2)
        
        if phase_retrieval_before(1) && do_phase_retrieval(1)
            fra_take_neg_log = 0;
        else
            fra_take_neg_log = 1;
        end
        fra_number_of_stds = 4;
        fra_vol_shape = [];
        
        if interactive_determination_of_rot_axis_slice > 1
            slice = interactive_determination_of_rot_axis_slice;
        elseif interactive_determination_of_rot_axis_slice <= 1 && interactive_determination_of_rot_axis_slice >= 0
            slice = round((raw_im_shape_binned2 - 1) * interactive_determination_of_rot_axis_slice + 1 );
        end
        fprintf( '\n slice : %u', slice)
        
        fprintf( '\n\nOFFSET:' )
        fprintf( '\n current rotation axis offset / position : %.2f, %.2f', rot_axis_offset, rot_axis_pos)
        fprintf( '\n calcul. rotation axis offset / position : %.2f, %.2f', rot_axis_offset_calc, rot_axis_pos_calc)
        fprintf( '\n default offset range : current ROT_AXIS_OFFSET + (-4:0.5:4)')
        offset = input( '\n\nENTER RANGE OF ROTATION AXIS OFFSETS (if empty use default range, if scalar skips interactive mode): ');
        if isempty( offset )
            % default range is centered at the given or calculated offset
            offset = rot_axis_offset + (-4:0.5:4);
        end
        
        % Loop over offsets
        while ~isscalar( offset )
            
            % Reco
            [vol, metrics_offset] = find_rot_axis_offset(proj, angles, slice, offset, rot_axis_tilt, fra_take_neg_log, fra_number_of_stds, fra_vol_shape, lamino, fixed_tilt, gpu_index, offset_shift);
            
            % Metric minima
            [~, min_pos] = min(cell2mat({metrics_offset(:).val}));
            [~, max_pos] = max(cell2mat({metrics_offset(:).val}));
            
            % Print image number, rotation axis values, and different metrics
            fprintf( '\n\nOFFSET:' )
            fprintf( '\n current rotation axis offset/position : %.2f, %.2f\n', rot_axis_offset, rot_axis_pos)
            fprintf( '%11s', 'image no.', 'offset', metrics_offset.name)
            for nn = 1:numel(offset)
                if offset(nn) == rot_axis_offset
                    cprintf( 'Green', sprintf('\n%11u%11.2f', nn, offset(nn)))
                else
                    cprintf( 'Black', '\n%11u%11.2f', nn, offset(nn))
                end
                
                for mm = 1:numel(metrics_offset)
                    if min_pos(mm) == nn
                        cprintf( 'Red', '%11.3g', metrics_offset(mm).val(nn) )
                    elseif max_pos(mm) == nn
                        cprintf( 'Blue', '%11.3g', metrics_offset(mm).val(nn) )
                    else
                        cprintf( 'Black', '%11.3g', metrics_offset(mm).val(nn) )
                    end
                end
            end
            
            % Plot metrics
            h_rot_off = figure('Name', 'OFFSET: metrics');
            x = [1:4 6:7];
            Y = cell2mat({metrics_offset(x).val});
            plot( offset, Y, '-+');
            axis tight
            legend( metrics_offset(x).name)
            title(sprintf('metric VS rotation axis offset'))
            drawnow
            
            % Play
            nimplay(vol, 1, [], 'OFFSET: sequence of reconstructed slices using different rotation axis offsets')
            
            % Input
            offset = input( '\n\nENTER ROTATION AXIS OFFSET OR A RANGE OF OFFSETS (if empty use current offset): ');
            if isempty( offset )
                offset = rot_axis_offset;
            end
            
            if isscalar( offset )
                fprintf( ' old rotation axis offset : %.2f', rot_axis_offset)
                rot_axis_offset = offset;
                fprintf( '\n new rotation axis offset : %.2f', rot_axis_offset)
                
                if interactive_determination_of_rot_axis_tilt(1)
                    
                    % Tilt
                    fprintf( '\n\nTILT:' )
                    fprintf( '\n current rotation axis tilt : %g rad = %g deg', rot_axis_tilt, rot_axis_tilt * 180 / pi)
                    fprintf( '\n calcul. rotation axis tilt : %g rad = %g deg', rot_axis_tilt_calc, rot_axis_tilt_calc * 180 / pi)
                    fprintf( '\n default tilt range is : current ROT_AXIS_TILT + (-0.005:0.001:0.005)')
                    tilt = input( '\n\nENTER TILT OF ROTATION AXIS OR RANGE OF TILTS (if empty use default):');
                    if isempty( tilt )
                        % default range is centered at the given or calculated tilt
                        tilt = rot_axis_tilt + (-0.005:0.001:0.005);
                    end
                    
                    while ~isscalar( tilt )
                        
                        % Reco
                        [vol, metrics_tilt] = find_rot_axis_tilt( proj, angles, slice, offset, tilt, fra_take_neg_log, fra_number_of_stds, fra_vol_shape, lamino, fixed_tilt, gpu_index, offset_shift);
                        
                        % Metric minima
                        [~, min_pos] = min(cell2mat({metrics_tilt(:).val}));
                        [~, max_pos] = max(cell2mat({metrics_tilt(:).val}));
                        
                        % Print image number and rotation axis tilt
                        fprintf( '%11s', 'image no.', 'tilt/rad', 'tilt/deg', metrics_tilt.name )
                        for nn = 1:numel(tilt)
                            if tilt(nn) == rot_axis_tilt
                                cprintf( 'Green', sprintf( '\n%11u%11g%11g', nn, tilt(nn), tilt(nn)/pi*180 ) )
                            else
                                cprintf( 'Black', sprintf( '\n%11u%11g%11g', nn, tilt(nn), tilt(nn)/pi*180 ) )
                            end
                            for mm = 1:numel(metrics_tilt)
                                if min_pos(mm) == nn
                                    cprintf( 'Red', '%11.3g', metrics_tilt(mm).val(nn) )
                                elseif max_pos(mm) == nn
                                    cprintf( 'Blue', '%11.3g', metrics_tilt(mm).val(nn) )
                                else
                                    cprintf( 'Black', '%11.3g', metrics_tilt(mm).val(nn) )
                                end
                            end
                        end
                        
                        % Plot metrics
                        h_rot_tilt = figure('Name', 'TILT: metrics');
                        x = 6:7;
                        Y = cell2mat({metrics_tilt(x).val});
                        plot( tilt, Y, '-+');
                        axis tight
                        legend( metrics_tilt(x).name)
                        title(sprintf('metric VS rotation axis tilt'))
                        drawnow
                        
                        % Play
                        nimplay(vol, 1, [], 'TILT: sequence of reconstructed slices using different rotation axis tilts')
                        
                        % Input
                        tilt = input( '\nENTER TILT OF ROTATION AXIS OR RANGE OF TILTS (if empty use current tilt):');
                        if isempty( tilt )
                            tilt = rot_axis_tilt;
                        end
                        
                        if isscalar( tilt )
                            rot_axis_tilt = tilt;
                            
                            rot_axis_pos = raw_im_shape_binned1 / 2 + rot_axis_offset;
                            
                            % Compare projection at 0 pi and projection at 1 pi corrected for rotation axis tilt
                            im1c = RotAxisSymmetricCropping( proj(:,rot_corr_area2,ind1), rot_axis_pos, 1);
                            im2c = flipud(RotAxisSymmetricCropping( proj(:,rot_corr_area2,ind2) , rot_axis_pos, 1));
                            [optimizer, metric] = imregconfig('monomodal');
                            tform_calc = imregtform(im2c, im1c, 'rigid', optimizer, metric);
                            rot_axis_tilt_calc = asin( tform_calc.T(1,2) ) / 2;
                            im2c_warped_calc =  imwarp(im2c, tform_calc, 'OutputView', imref2d(size(im1c)));
                            
                            tform_cur = tform_calc;
                            %tform_cur.Dimensionality = 2;
                            tform_cur.T = [cos( 2 * rot_axis_tilt ) sin( 2 * rot_axis_tilt ) 0; ...
                                -sin( 2 * rot_axis_tilt ) cos( 2 * rot_axis_tilt ) 0 ; ...
                                tform_calc.T(3,1) tform_calc.T(3,2) 1];
                            %0 0 1];
                            
                            im2c_warped_cur =  imwarp(im2c, tform_cur, 'OutputView', imref2d(size(im1c)));
                            
                            x = ceil( 2 * abs( sin(rot_axis_tilt) ) * max( size(im1c)) ) + 2;
                            
                            im2c_warped_cur(im2c_warped_cur == 0) = mean2( im2c );
                            im2c_warped_calc(im2c_warped_calc == 0) = mean2( im2c );
                            
                            fprintf( '\n current rotation axis tilt: %g rad (%g deg)', rot_axis_tilt, rot_axis_tilt * 180 / pi)
                            fprintf( '\n calcul. rotation axis tilt: %g rad (%g deg)', rot_axis_tilt_calc, rot_axis_tilt_calc * 180 / pi)
                            
                            name = sprintf( 'projections at 0 and 180 degree. corrected. CURRENT rot axis tilt: %g, rot axis offset: %g', rot_axis_tilt, rot_axis_offset);
                            nimplay( cat(3, im1c(x:end-x,x:end-x)', im2c_warped_cur(x:end-x,x:end-x)'), 1, 0, name)
                            
                            name = sprintf( 'projections at 0 and 180 degree. corrected. CALCULATED rot axis tilt: %g, rot axis offset: %g', rot_axis_tilt_calc, rot_axis_offset);
                            nimplay( cat(3, im1c(x:end-x,x:end-x)', im2c_warped_calc(x:end-x,x:end-x)'), 1, 0, name)
                            
                            tilt = input( '\nENTER ROTATION AXIS TILT (if empty use current value): ');
                            if isempty( tilt )
                                tilt = rot_axis_tilt;
                            else
                                rot_axis_tilt = tilt;
                            end
                            
                            offset = input( '\nENTER RANGE OF OFFSETS TO CONTINUE INTERACTIVE LOOP OR TYPE ENTER TO EXIT LOOP: ');
                            if isempty( offset )
                                offset = rot_axis_offset;
                            end
                        end
                    end
                else
                    rot_axis_tilt = 0;
                end
            end
        end
        
        rot_axis_pos = raw_im_shape_binned1 / 2 + rot_axis_offset;
        fprintf( '\nEND OF INTERACTIVE MODE' )
        
        tint = toc - tint;
    end
    PrintVerbose(verbose, '\n rotation axis offset: %.2f', rot_axis_offset );
    PrintVerbose(verbose, '\n rotation axis position: %.2f', rot_axis_pos );
    PrintVerbose(verbose, '\n rotation axis tilt: %g rad (%g deg)', rot_axis_tilt, rot_axis_tilt * 180 / pi)
    
    if interactive_determination_of_rot_axis_tilt(1) && visual_output(1)
        h4 = figure('Name','Projections at 0 and pi cropped symmetrically to rotation center');
        n = 2;
        m = 2;
        
        subplot(m, n, 1)
        imsc1( im1c )
        axis equal tight
        title(sprintf('proj at 0'))
        colorbar
        
        subplot(m, n, 2)
        imsc1( im2c )
        axis equal tight
        title(sprintf('proj at pi'))
        colorbar
        
        x = ceil(1 + 2 * abs( sin(rot_axis_tilt) ) * max( size(im1c)) );
        
        subplot(m, n, 3)
        imsc1( abs( im1c(x:end-x,x:end-x) - im2c(x:end-x,x:end-x) ) )
        axis equal tight
        title(sprintf('difference before'))
        colorbar
        
        subplot(m, n, 4)
        
        %% tform_cal or tform_cur
        [optimizer, metric] = imregconfig('monomodal');
        tform_calc = imregtform(im2c, im1c, 'rigid', optimizer, metric);
        im2c_warped_cur =  imwarp(im2c, tform_calc, 'OutputView',imref2d(size(im1c)));
        imsc1( abs( im1c(x:end-x,x:end-x) - im2c_warped_cur(x:end-x,x:end-x) ) )
        axis equal tight
        title(sprintf('difference corrected'))
        colorbar
        
        drawnow
    end
end

%% Stitch projections
if stitch_projections(1)
    t = toc;
    PrintVerbose(verbose, '\nStitch projections:')
    if rot_angle_full < 1.9 * pi
        error( 'full angle of rotation smaller than 2 pi: %g pi', rot_angle_full/pi)
    end
    % number of unstitched projections
    num_proj = size( proj, 3);
    % last projection within [0,pi)
    [~, num_proj_sti] = min( abs(angles - pi));
    % number of stitched projections
    num_proj_sti = num_proj_sti - 1;
    % index range of projections to be stitched
    xl = 1:round(rot_axis_pos);
    xr = 1:xl(end)-1;
    im_shape_sti1 = numel( xl ) + numel( xr );
    % Preallocation
    proj_sti = zeros( im_shape_sti1 , raw_im_shape_binned2, num_proj_sti, 'single');
    for nn = 1:num_proj_sti
        nn2 = mod(num_proj_sti + nn - 1, num_proj) + 1;
        im = zeros( im_shape_sti1, raw_im_shape_binned2);
        switch lower( stitch_method )
            case 'step'
                im = cat(1, proj(xl,:,nn), flipud( proj(xr,:,nn2) ) );
            case {'linear', 'sine'}
                % overlap region
                overlap = 2 * rot_axis_pos - raw_im_shape_binned1 : raw_im_shape_binned1;
                % overlap ramp
                x = (0:1/(numel(overlap)-1):1);
                % 1D weight
                w = ones(raw_im_shape_binned1, 1);
                switch lower( stitch_method )
                    case 'linear'
                        w(overlap) = 1 - x;
                    case 'sine'
                        w(overlap) = 0.5 * cos(pi*x) + 0.5;
                end
                % weighted projections
                iml = bsxfun(@times, proj(:,:,nn), w);
                imr = flipud( bsxfun(@times, proj(:,:,nn2), w ) );
                % stitched projection
                im(1:raw_im_shape_binned1,:) = iml;
                im(end - raw_im_shape_binned1 + 1:end,:) = im(end - raw_im_shape_binned1 + 1:end,:) + imr;
        end
        proj_sti(:,:,nn) = im;
    end
    pause(0.01)
    proj = proj_sti;
    clear proj_sti;
    angles = angles(1:num_proj_sti);
    PrintVerbose(verbose, ' Time elapsed: %.1f (%.2f min)', toc-t, (toc-t)/60)
    PrintVerbose(verbose, '\n shape of stitched projections : %u %u %u', size( proj ) )
    PrintVerbose(verbose, '\n memory allocated : %.2f GiB', Bytes( proj, 3 ) )
end

%% Crop projections at rotation axis position
% This is to avoid oversampling for scans with excentric rotation axis and
% reconstructing WITHOUT stitching
if crop_at_rot_axis(1)
    switch excentric_rot_axis
        case 1
            proj( ceil(rot_axis_pos) + 1:end, :, :) = [];
        case -1
            %% CHECK
            proj( 1:floor(rot_axis_pos)-1, :, :) = [];
    end
    if isempty( vol_shape )
        vol_shape = [raw_im_shape_binned1, raw_im_shape_binned1, raw_im_shape_binned2];
    end
end

%% Save sinograms
if write_sino(1)
    t = toc;
    PrintVerbose(verbose, '\nSave sinogram:')
    CheckAndMakePath(sino_path)
    % Save slices
    parfor nn = 1:raw_im_shape_binned2
        filename = sprintf( '%ssino_%06u.tif', sino_path, nn);
        write32bitTIFfromSingle( filename, squeeze( proj( :, nn, :) )' )
    end
    pause(0.01)
    PrintVerbose(verbose, ' done.')
    PrintVerbose(verbose, ' Time elapsed: %.1f s (%.2f min)', toc-t, (toc-t)/60)
end

%% Phase retrieval %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
edp = [energy, sample_detector_distance, eff_pixel_size_binned];
if do_phase_retrieval(1) && ~phase_retrieval_before(1)
    PrintVerbose(verbose, '\nPhase retrieval')
    t = toc;
    PrintVerbose(verbose, '\n energy : %g eV', energy)
    PrintVerbose(verbose, '\n sample detector distance : %g m', sample_detector_distance)
    
    if isempty( take_neg_log )
        take_neg_log = 0;
    end
    
    % Phase retrieval filter
    im_shape = [size(proj,1) , size(proj,2)];
    im_shape_pad = (1 + phase_padding) * im_shape;
    [phase_filter, pha_appendix] = PhaseFilter( phase_retrieval_method, im_shape_pad, edp, phase_retrieval_reg_par, phase_retrieval_bin_filt, phase_retrieval_cutoff_frequ, 'single');
    
    % reco phase dir
    if isempty( subfolder_reco )
        reco_phase_path = [out_path, filesep, 'reco_phase', filesep, pha_appendix, filesep];
    else
        reco_phase_path = [out_path, filesep, 'reco_phase', filesep, pha_appendix, filesep, subfolder_reco, filesep];
    end
    CheckAndMakePath( reco_phase_path )
    PrintVerbose(verbose, '\n reco_phase_path : %s', reco_phase_path)
    PrintVerbose(verbose, '\n phase retrieval method : %s', phase_retrieval_method)
    
    % Retrieval
    parfor nn = 1:size(proj, 3)
        % combined GPU and parfor usage requires memory management
        im = padarray( proj(:,:,nn), phase_padding * im_shape, 'symmetric', 'post' );
        %im = padarray( gpuArray( proj(:,:,nn) ), raw_im_shape_binned, 'post', 'symmetric' );
        im = -real( ifft2( phase_filter .* fft2( im ) ) );
        proj(:,:,nn) = im(1:im_shape(1), 1:im_shape(2));
        %proj(:,:,nn) = gather( im(1:raw_im_shape_binned1, 1:raw_im_shape_binned2) );
    end
    pause(0.01)
    PrintVerbose(verbose, '\n done in %g s (%.2f min)', toc-t, (toc-t)/60)
    
    % Save phase maps
    if write_phase_map(1)
        t = toc;
        PrintVerbose(verbose, '\nSave phase maps:')
        phase_map_path = [phase_map_path, pha_appendix, filesep];
        CheckAndMakePath( phase_map_path );
        parfor nn = 1:size( proj, 3)
            filename = sprintf( '%sphase_%06u.tif', phase_map_path, nn);
            write32bitTIFfromSingle( filename, squeeze( rot90( proj( :, :, nn) ) ) )
        end
        pause(0.01)
        PrintVerbose(verbose, ' Time elapsed: %.1f s (%.2f min)', toc-t, (toc-t)/60)
    end
end

%% Save sinograms of phase maps
if do_phase_retrieval(1) && write_sino_phase(1)
    t = toc;
    PrintVerbose(verbose, '\nSave phase map sinogram:')
    CheckAndMakePath(sino_phase_path)
    % Save slices
    parfor nn = 1:raw_im_shape_binned2
        filename = sprintf( '%ssino_%06u.tif', sino_phase_path, nn);
        write32bitTIFfromSingle( filename, squeeze( proj( :, nn, :) )' )
    end
    pause(0.01)
    PrintVerbose(verbose, ' done.')
    PrintVerbose(verbose, ' Time elapsed: %.1f s (%.2f min)', toc - t, (toc - t) / 60)
end

%% Post phase retrieval binning
if do_phase_retrieval(1) && ( phase_bin(1) > 1)
    t = toc;
    PrintVerbose(verbose, 'Post phase retrieval binning:')
    proj_bin = zeros( floor(size( proj ) ./ [2 2 1] ), 'single');
    parfor nn = 1:size( proj, 3)
        proj_bin(:,:,nn) = Binning( proj(:,:,nn), phase_bin ) / phase_bin^2;
    end
    proj = proj_bin;
    clear proj_bin;
    rot_axis_pos = rot_axis_pos / phase_bin;
    rot_axis_offset = rot_axis_offset / phase_bin;
    PrintVerbose(verbose, ' done in %g s (%.2f min)', toc - t, (toc - t) / 60)
end

%% Tomographic reco
if do_tomo(1)
    PrintVerbose(verbose, '\nTomographic reconstruction:')
    
    if stitch_projections(1)
        rot_axis_offset_reco = 0;
    elseif crop_at_rot_axis(1)
        rot_axis_offset_reco = rot_axis_pos - size( proj, 1) / 2;
    else
        rot_axis_offset_reco = rot_axis_offset;
    end
    
    %% Volume shape and size
    if isempty( vol_shape )
        % default volume given by the detector width and height
        vol_shape1 = size( proj, 1);
        vol_shape3 = size( proj, 2 );
        vol_shape = [vol_shape1, vol_shape1, vol_shape3];
    else
        vol_shape1 = size( proj, 1);
        if vol_shape(1) <=  1
            vol_shape(1) = round( 1 + vol_shape(1) * (vol_shape1 - 1) );
        end
        if vol_shape(2) <=  1
            vol_shape(2) = round( 1 + vol_shape(2) * (vol_shape1 - 1) );
        end
        if vol_shape(3) <=  1
            vol_shape(3) = round( 1 + vol_shape(3) * (vol_shape3 - 1) );
        end
    end
    if isempty( vol_size )
        vol_size = [-vol_shape(1)/2, vol_shape(1)/2, -vol_shape(2)/2, vol_shape(2)/2, -vol_shape(3)/2, vol_shape(3)/2];
    end
    PrintVerbose(verbose, '\n shape of reconstructed volume: [%g, %g, %g]', vol_shape )
    PrintVerbose(verbose, '\n memory required: %.2f GiB', prod( vol_shape ) * 4 / 1024^3 )
    
    if isempty( take_neg_log )
        take_neg_log = 1;
    end
    angles_reco = angles;
    % Delete redundant projection
    if isequal( angles(1), angles(end) )
        angles_reco(end) = [];
        proj(:,:,end) = [];
    end
    
    % Filter sinogram
    PrintVerbose(verbose, '\n Filter sino:' )
    t2 = toc;
    filt = iradonDesignFilter(fbp_filter_type, (1 + fbp_filter_padding) * size( proj, 1), fpb_filter_freq_cutoff);
    if butterworth_filter(1)
        [b, a] = butter(butterworth_order, butterworth_cutoff_frequ);
        bw = freqz(b, a, numel(filt) );
        filt = filt .* bw;
    end
    
    proj_shape1 = size( proj, 1);
    parfor nn =  1:size( proj, 2)
        im = proj(:,nn,:);
        im = padarray( NegLog(im, take_neg_log), fbp_filter_padding * [proj_shape1 0 0], fbp_filter_padding_method, 'post' );
        im = real( ifft( bsxfun(@times, fft( im, [], 1), filt), [], 1, 'symmetric') );
        proj(:,nn,:) = im(1:proj_shape1,:,:);
    end
    pause(0.01)
    PrintVerbose(verbose, ' done in %.2f min.', (toc - t2) / 60)
    
    if crop_at_rot_axis(1)
        % half weight pixel at rot axis pos as it is used twice
        switch excentric_rot_axis
            case 1
                proj( end, :, :) = 0.5 * proj( end, :, :) ;
            case -1
                proj( 1, :, :) = 0.5 * proj( 1, :, :) ;
        end
    end
    
    % Backprojection
    PrintVerbose(verbose, '\n Backproject:')
    t2 = toc;
    vol = astra_parallel3D( permute(proj, [1 3 2]), rot_angle_offset + angles_reco, rot_axis_offset_reco + offset_shift, vol_shape, vol_size, astra_pixel_size, link_data, rot_axis_tilt, gpu_index);
    pause(0.01)
    PrintVerbose(verbose, ' done in %.2f min.', (toc - t2) / 60)
    
    vol_min = min( vol(:) );
    vol_max = max( vol(:) );
    
    % Show orthogonal vol cuts
    if visual_output(1)
        figure('Name', 'Volum cuts');
        
        subplot(1,3,1)
        nn = round( size( vol, 3 ) / 2);
        imsc( squeeze( vol(:,:,nn) ) )
        axis equal tight
        title( sprintf( 'vol z = %u', nn ) )
        colorbar
        
        subplot(1,3,2)
        nn = round( size( vol, 2 ) / 2);
        imsc( squeeze( vol(:,nn,:) ) )
        axis equal tight
        title( sprintf( 'vol y = %u', nn ) )
        colorbar
        
        subplot(1,3,3)
        nn = round( size( vol, 1 ) / 2);
        imsc( squeeze( vol(nn,:,:) ) )
        axis equal tight
        title( sprintf( 'vol x = %u', nn ) )
        colorbar
        
        drawnow
    end
    
    % Save volume
    if write_reco(1)
        if do_phase_retrieval(1)
            reco_path = reco_phase_path;
        end
        CheckAndMakePath( reco_path )
        
        % Save reco path to file
        filename = [userpath, filesep, 'experiments/p05/path_to_latest_reco'];
        fid = fopen( filename , 'w' );
        fprintf( fid, '%s', reco_path );
        fclose( fid );
        
        % Single precision: 32-bit float tiff
        write_volume( write_float, vol, 'float', reco_path, raw_bin, phase_bin, 1, 0, verbose);
        
        % Compression of dynamic range
        if write_8bit || write_8bit_binned || write_16bit || write_16bit_binned
            [tlow, thigh] = compression( vol, compression_method, compression_parameter );
        else
            tlow = 0;
            thigh = 1;
        end
        
        % 16-bit tiff
        write_volume( write_16bit, (vol - tlow)/(thigh - tlow), 'uint16', reco_path, raw_bin, phase_bin, 1, 0, verbose);
        
        % 8-bit tiff
        write_volume( write_8bit, (vol - tlow)/(thigh - tlow), 'uint8', reco_path, raw_bin, phase_bin, 1, 0, verbose);
        
        % Bin data
        if write_float_binned || write_16bit_binned || write_8bit_binned || write_8bit_segmented
            PrintVerbose(verbose, '\n Binning:')
            t2 = toc;
            vol = Binning( vol, reco_bin ) / reco_bin^3;
            PrintVerbose(verbose, ' done in %.2f min.', (toc - t2) / 60)
        end
        
        % Binned single precision: 32-bit float tiff
        write_volume( write_float_binned, vol, 'float', reco_path, raw_bin, phase_bin, reco_bin, 0, verbose);
        
        % 16-bit tiff binned
        write_volume( write_16bit_binned, (vol - tlow)/(thigh - tlow), 'uint16', reco_path, raw_bin, phase_bin, reco_bin, 0, verbose);
        
        % 8-bit tiff binned
        write_volume( write_8bit_binned, (vol - tlow)/(thigh - tlow), 'uint8', reco_path, raw_bin, phase_bin, reco_bin, 0, verbose);
        
        % segmentation
        if write_8bit_segmented(1)
            [vol, out] = segment_volume(vol, 2^10, visual_output, verbose);
            save_path = write_volume( 1, vol/255, 'uint8', reco_path, raw_bin, phase_bin, reco_bin, 0, verbose, '_segmented');
            save( sprintf( '%ssegmentation_info.m', save_path), 'out', '-mat', '-v7.3')
        end
        
    end
    
    PrintVerbose(verbose, ' Time elapsed: %.1f s (%.2f min)', toc-t, (toc-t)/60 )
end

%% Log file
if write_reco(1)
    logfile_path = reco_path;
else
    logfile_path = out_path;
end
logfile_name = sprintf( '%sreco.log', logfile_path);
fid = fopen( logfile_name, 'w');
fprintf(fid, 'scan_name : %s\n', scan_name);
fprintf(fid, 'beamtime_id : %s\n', beamtime_id);
fprintf(fid, 'scan_path : %s\n', scan_path);
fprintf(fid, 'reco_path : %s\n', reco_path);
fprintf(fid, 'MATLAB notation, index of first element: 1, range: first:stride:last\n');
fprintf(fid, 'MATLAB version : %s\n', version);
fprintf(fid, 'platform : %s\n', computer);
fprintf(fid, 'camera : %s\n', cam);
fprintf(fid, 'num_dark_found : %u\n', num_dark);
fprintf(fid, 'num_ref_found : %u\n', num_ref_found);
fprintf(fid, 'num_ref_used : %u\n', num_ref_used);
fprintf(fid, 'ref_range : %u:%u:%u\n', ref_range(1), ref_range(2) - ref_range(1), ref_range(end) );
fprintf(fid, 'num_proj_found : %u\n', num_proj_found);
fprintf(fid, 'num_proj_used : %u\n', num_proj_used);
fprintf(fid, 'proj_range : %u:%u:%u\n', proj_range(1), proj_range(2) - proj_range(1), proj_range(end) );
fprintf(fid, 'raw_im_shape_uncropped : %u %u\n', raw_im_shape_uncropped);
fprintf(fid, 'raw_roi : %u %u\n', raw_roi);
fprintf(fid, 'raw_image_shape : %u %u\n', raw_im_shape);
fprintf(fid, 'raw_image_shape_binned : %u %u\n', raw_im_shape_binned);
fprintf(fid, 'raw_binning_factor : %u\n', raw_bin);
fprintf(fid, 'bin_before_filtering : %u\n', bin_before_filtering);
fprintf(fid, 'effective_pixel_size : %g micron\n', eff_pixel_size * 1e6);
fprintf(fid, 'effective_pixel_size_binned : %g micron\n', eff_pixel_size_binned * 1e6);
fprintf(fid, 'energy : %g eV\n', energy);
fprintf(fid, 'sample_detector_distance : %f m\n', sample_detector_distance);
fprintf(fid, 'ring_current_normalization : %u\n', ring_current_normalization);
fprintf(fid, 'proj_flat_correlation_method : %s\n', correlation_method);
fprintf(fid, 'proj_flat_correlation_num_flats : %u\n', corr_num_flats);
fprintf(fid, 'flat_field_correlation_area_1 : %u:%u:%u\n', flat_corr_area1(1), flat_corr_area1(2) - flat_corr_area1(1), flat_corr_area1(end));
fprintf(fid, 'flat_field_correlation_area_2 : %u:%u:%u\n', flat_corr_area2(1), flat_corr_area2(2) - flat_corr_area2(1), flat_corr_area2(end));
if ~read_flatcor
    fprintf(fid, 'min_max_of_all_darks : %6g %6g\n', darks_min, darks_max);
    fprintf(fid, 'min_max_of_median_dark : %6g %6g\n', dark_med_min, dark_med_max);
    fprintf(fid, 'min_max_of_all_flats : %6g %6g\n', flat_min, flat_max);
    fprintf(fid, 'min_max_of_all_corrected_flats : %6g %6g\n', flat_min2, flat_max2);
    fprintf(fid, 'min_max_of_all_raws :  %6g %6g\n', raw_min, raw_max);
    fprintf(fid, 'min_max_of_all_corrected_raws :  %6g %6g\n', raw_min2, raw_max2);
    fprintf(fid, 'min_max_of_all_flat_corr_projs : %g %g \n', proj_min, proj_max);
end
% Phase retrieval
if do_phase_retrieval(1)
    fprintf(fid, 'do_phase_retrieval : %u\n', do_phase_retrieval);
    fprintf(fid, 'phase_bin : %u\n', phase_bin);
    fprintf(fid, 'phase_retrieval_method : %s\n', phase_retrieval_method);
    fprintf(fid, 'phase_retrieval_regularisation_parameter : %f\n', phase_retrieval_reg_par);
    fprintf(fid, 'phase_retrieval_binary_filter_threshold : %f\n', phase_retrieval_bin_filt);
    fprintf(fid, 'phase_retrieval_cutoff_frequ : %f pi\n', phase_retrieval_cutoff_frequ / pi);
    fprintf(fid, 'phase_padding : %u\n', phase_padding);
end
% Volume
fprintf(fid, 'volume_shape : %u %u %u\n', vol_shape(1), vol_shape(2), vol_shape(3));
fprintf(fid, 'volume_size : %f %f %f %f %f %f\n', vol_size(1), vol_size(2), vol_size(3), vol_size(4), vol_size(5), vol_size(6));
% Rotation
fprintf(fid, 'excentric_rot_axis : %i\n', excentric_rot_axis);
fprintf(fid, 'crop_at_rot_axis : %u\n', crop_at_rot_axis);
fprintf(fid, 'stitch_projections : %u\n', stitch_projections);
fprintf(fid, 'stitch_method : %s\n', stitch_method );
fprintf(fid, 'rotation_angle_full_rad : %f * pi\n', rot_angle_full / pi);
fprintf(fid, 'rotation_angle_offset_rad : %f * pi\n', rot_angle_offset / pi);
fprintf(fid, 'rotation_axis_offset_calculated : %f\n', rot_axis_offset_calc);
fprintf(fid, 'rotation_axis_offset : %f\n', rot_axis_offset);
fprintf(fid, 'rotation_axis_offset_reco : %f\n', rot_axis_offset_reco);
fprintf(fid, 'rotation_axis_position_calculated : %f\n', rot_axis_pos_calc);
fprintf(fid, 'rotation_axis_position : %f\n', rot_axis_pos);
fprintf(fid, 'rotation_axis_tilt_calculated : %f\n', rot_axis_tilt_calc);
fprintf(fid, 'rotation_axis_tilt : %f\n', rot_axis_tilt);
fprintf(fid, 'raw_image_binned_center : %f\n', raw_im_shape_binned1 / 2);
fprintf(fid, 'rotation_correlation_area_1 : %u:%u:%u\n', rot_corr_area1(1), rot_corr_area1(2) - rot_corr_area1(1), rot_corr_area1(end));
fprintf(fid, 'rotation_correlation_area_2 : %u:%u:%u\n', rot_corr_area2(1), rot_corr_area2(2) - rot_corr_area2(1), rot_corr_area2(end));
fprintf(fid, 'interactive_determination_of_rot_axis : %u\n', interactive_determination_of_rot_axis);
% FBP
fprintf(fid, 'ring_filter : %u\n', ring_filter);
fprintf(fid, 'ring_filter_method : %s\n', ring_filter_method);
fprintf(fid, 'ring_filter_jm_median_width : %s\n', sprintf( '%u ', ring_filter_jm_median_width) );
fprintf(fid, 'fbp_filter_type : %s\n', fbp_filter_type);
fprintf(fid, 'fpb_filter_freq_cutoff : %f\n', fpb_filter_freq_cutoff);
fprintf(fid, 'fbp_filter_padding : %u\n', fbp_filter_padding);
fprintf(fid, 'fbp_filter_padding_method : %s\n', fbp_filter_padding_method);
fprintf(fid, 'butterworth_filter : %u\n', butterworth_filter);
fprintf(fid, 'butterworth_order : %u\n', butterworth_order);
fprintf(fid, 'butterworth_cutoff_frequency : %f\n', butterworth_cutoff_frequ);
fprintf(fid, 'astra_pixel_size : %f\n', astra_pixel_size);
fprintf(fid, 'take_negative_logarithm : %u\n', take_neg_log);
fprintf(fid, 'gpu_name : %s\n', gpu.Name);
if exist( 'vol_min', 'var')
    fprintf(fid, 'volume_min-max : %g %g\n', vol_min, vol_max);
end
fprintf(fid, 'write_float : %u\n', write_float);
fprintf(fid, 'write_16bit : %g\n', write_16bit);
fprintf(fid, 'write_8bit : %u\n', write_8bit);
fprintf(fid, 'write_8bit_binned : %u\n', write_8bit_binned);
fprintf(fid, 'compression_method : %s\n', compression_method);
fprintf(fid, 'compression_limits : %f %f\n', tlow, thigh);
fprintf(fid, 'reco_bin : %u\n', reco_bin);
fprintf(fid, 'full_reconstruction_time : %.1f s\n', toc);
fprintf(fid, 'date_of_reconstruction : %s\n', datetime);
fprintf(fid, 'rot_axis_offset at %u x binning : %f\n', raw_bin, rot_axis_offset);
fclose(fid);
% End of log file

PrintVerbose(verbose, '\n log file : %s', logfile_name)
PrintVerbose(verbose, '\n reco_path : \n%s', reco_path)
PrintVerbose(verbose && interactive_determination_of_rot_axis, '\nTime elapsed in interactive mode: %g s (%.2f min)', tint, tint / 60 );
PrintVerbose(verbose, '\nTime elapsed for computation: %g s (%.2f min)', toc - tint, (toc - tint) / 60 );
PrintVerbose(verbose, '\nFINISHED: %s\n\n', scan_name)
% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dbclear if error
