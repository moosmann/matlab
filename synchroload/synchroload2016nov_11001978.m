function synchroload2016nov_11001978(SUBSETS, RUN_RECO, PRINT_PARAMETERS )
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

read_flatcor = 0; % read flatfield-corrected images from disc, skips preprocessing
read_flatcor_path = ''; % subfolder of 'flat_corrected' containing projections
% PREPROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_roi = [];-1; % if []: no ROI; if [y0 y1]: vertical ROI, skips first raw_roi(1)-1 lines, reads until raw_roi(2). When raw_roi(2) < 0 reads until end - |raw_roi(2)|; if negative scalar: selects ROI automatically.Not working for *.raw data where images are flipped.
raw_bin = 2; % projection binning factor: 1, 2, or 4
bin_before_filtering = 0; % Binning is applied before pixel filtering, much faster but less effective filtering
excentric_rot_axis = 0; % off-centered rotation axis increasing FOV. -1: left, 0: centeerd, 1: right. influences rot_corr_area1
crop_at_rot_axis = 0; % for recos of scans with excentric rotation axis but WITHOUT projection stitching
stitch_projections = 0; % for 2 pi scans: stitch projection at rotation axis position. Recommended with phase retrieval to reduce artefacts. Standard absorption contrast data should work well without stitching. Subpixel stitching not supported (non-integer rotation axis position is rounded, less/no binning before reconstruction can be used to improve precision).
stitch_method = 'sine'; 'step';'linear'; %  ! CHECK correlation area !
proj_range = 1; % range of projections to be used (from all found, if empty or 1: all, if scalar: stride, if range: start:incr:end
ref_range = 1; % range of flat fields to be used (from all found), if empty or 1: all. if scalar: stride, if range: start:incr:end
energy = []; % in eV! if empty: read from log file
sample_detector_distance = []; % in m. if empty: read from log file
eff_pixel_size = []; % in m. if empty: read from log file. effective pixel size =  detector pixel size / magnification
correlation_method = 'ssim-ml';'none';'entropy';'diff';'shift';'ssim';'std';'cov';'corr';'cross-entropy12';'cross-entropy21';'cross-entropyx';'none';
corr_shift_max_pixelshift = 0.25; % maximum pixelshift allowed for 'shift'-correlation method: if 0 use the best match (i.e. the one with the least shift), if > 0 uses all flats with shifts smaller than corr_shift_max_pixelshift
corr_num_flats = 1; % number of flat fields used for average/median of flats. for 'shift'-correlation its the maximum number
ring_current_normalization = 1; % normalize flat fields and projections by ring current
flat_corr_area1 = [0 0.02];%[0.98 1];% correlation area: index vector or relative/absolute position of [first pix, last pix]
flat_corr_area2 = [0.2 0.8]; % correlation area: index vector or relative/absolute position of [first pix, last pix]
decimal_round_precision = 2; % precision when rounding pixel shifts
strong_abs_thresh = 1; % if 1: does nothing, if < 1: flat-corrected values below threshold are set to one
ring_filter.apply = 0; % ring artifact filter (use only for scans without lateral sample movement)
ring_filter.apply_before_stitching = 0; % ! Consider when phase retrieval is applied !
ring_filter.method = 'wavelet-fft';'jm'; 
ring_filter.waveletfft.dec_levels = 2:5; % decomposition levels for 'wavelet-fft'
ring_filter.waveletfft.wname = 'db25';'db30'; % wavelet type for 'wavelet-fft'
ring_filter.waveletfft.sigma = 2.4; %  suppression factor for 'wavelet-fft'
ring_filter.jm.median_width = 11; % multiple widths are applied consecutively, eg [3 11 21 31 39];
% PHASE RETRIEVAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phase_retrieval.apply = 0; % See 'PhaseFilter' for detailed description of parameters !
phase_retrieval.apply_before = 0; % before stitching, interactive mode, etc. For phase-contrast data with an excentric rotation axis phase retrieval should be done afterwards. To find the rotataion axis position use this option in a first run, and then turn it of afterwards.
phase_retrieval.post_binning_factor = 1; % Binning factor after phase retrieval, but before tomographic reconstruction
phase_retrieval.method = 'tie';'qpcut'; %'qp' 'ctf' 'tie' 'qp2' 'qpcut'
phase_retrieval.reg_par = 1.5; % regularization parameter. larger values tend to blurrier images. smaller values tend to original data.
phase_retrieval.bin_filt = 0.15; % threshold for quasiparticle retrieval 'qp', 'qp2'
phase_retrieval.cutoff_frequ = 2 * pi; % in radian. frequency cutoff in Fourier space for 'qpcut' phase retrieval
phase_retrieval.padding = 1; % padding of intensities before phase retrieval, 0: no padding
% TOMOGRAPHY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
do_tomo = 1; % run tomographic reconstruction
vol_size = [];%[-0.6 0.6 -0.6 0.6 -0.5 0.5];% for excentric rot axis pos; 6-component vector [xmin xmax ymin ymax zmin zmax]. if empty, volume is centerd within vol_shape. unit voxel size is assumed. if smaller than 10 values are interpreted as relative size w.r.t. the detector size. Take care bout minus signs!
vol_shape = [];%[1.2 1.2 1];% for excentric rot axis pos; % shape (voxels) of reconstruction volume. in absolute number of voxels or in relative number w.r.t. the default volume which is given by the detector width and height
rot_angle_full = []; % in radians: empty ([]), full angle of rotation, or array of angles. if empty full rotation angles is determined automatically to pi or 2 pi
rot_angle_offset = 1*pi; % global rotation of reconstructed volume
rot_axis_offset = 0; % if empty automatic computation is used (not very reliable)
rot_axis_pos = []; % if empty use automatic computation. EITHER OFFSET OR POSITION MUST BE EMPTY. YOU MUST NOT USE BOTH!
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
butterworth_cutoff_frequ = 0.9;
astra_pixel_size = 1; % detector pixel size for reconstruction: if different from one 'vol_size' must to be ajusted, too!
take_neg_log = []; % take negative logarithm. if empty, use 1 for attenuation contrast, 0 for phase contrast
% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out_path = '/gpfs/petra3/scratch/moosmanj';% absolute path were output data will be stored. !!overwrites the write.to_scratch flag. if empty uses the beamtime directory and either 'processed' or 'scratch_cc'
write.to_scratch = 0; % write to 'scratch_cc' instead of 'processed'
write.flatcor = 1; % save preprocessed flat corrected projections
write.phase_map = 0; % save phase maps (if phase retrieval is not 0)
write.sino = 0; % save sinograms (after preprocessing & before FBP filtering and phase retrieval)
write.phase_sino = 0; % save sinograms of phase maps
write.reco = 1; % save reconstructed slices (if do_tomo=1)
write.float = 1; % single precision (32-bit float) tiff
write.uint16 = 0;
write.uint8 = 0;
write.post_reconstruction_binning_factor = 2; % IF BINNED VOLUMES ARE SAVED: binning factor of reconstructed volume
write.float_binned = 1; % binned single precision (32-bit float) tiff
write.uint16_binned = 0;
write.uint8_binned = 0;
parfolder = '';% parent folder for 'reco', 'sino', 'phase', and 'flat_corrected'
subfolder_flatcor = ''; % subfolder in 'flat_corrected'
subfolder_phase_map = ''; % subfolder in 'phase_map'
subfolder_sino = ''; % subfolder in 'sino'
subfolder_reco = ''; % subfolder in 'reco'
write.uint8_segmented = 0;
% INTERACTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
verbose = 1; % print information to standard output
visual_output = 1; % show images and plots during reconstruction
interactive_determination_of_rot_axis = 0; % reconstruct slices with dif+ferent rotation axis offsets
interactive_determination_of_rot_axis_tilt = 0; % reconstruct slices with different offset AND tilts of the rotation axis
SET_DEFAULT

%% Data sets %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_path = '/asap3/petra3/gpfs/p05/2016/data/11001978/raw/';

% corroded screw
scan_path = [raw_path 'mah_01'];
rot_axis_offset = -135.75 / raw_bin;
rot_axis_tilt = -0.003;
%write_sino = 1;
ADD

% corroded screw
scan_path = [raw_path 'mah_02'];
rot_axis_offset = -135.75 / raw_bin;
rot_axis_tilt = -0.003;
ADD

% corroded screw
scan_path = [raw_path 'mah_03'];
rot_axis_offset = -135.75 / raw_bin;
rot_axis_tilt = -0.003;
ADD

% implant fresh
scan_path = [raw_path 'mah_04'];
excentric_rot_axis = 1;
crop_at_rot_axis = 1;
rot_axis_offset = 628 / raw_bin;
rot_axis_tilt = -0.003;
%vol_shape = [2155 2155 1050];
ADD('r')

% corroded screw
scan_path = [raw_path 'mah_05'];
rot_axis_offset = 2 / raw_bin;
rot_axis_tilt = -0.003;
ADD

% corroded screw
scan_path = [raw_path 'mah_06_Mg10G_004'];
rot_axis_offset = 5 / raw_bin;
rot_axis_tilt = -0.003;
ADD

% implant fresh formalin
scan_path = [raw_path 'mah_07_bone_in_formalin'];
rot_axis_offset = 2 / raw_bin;
rot_axis_tilt = -0.003;
ADD

% corrosion cell
scan_path = [raw_path 'mah_08_corrosion_cell_A'];
rot_axis_offset = -3.5 / raw_bin;
rot_axis_tilt = -0.003;
ADD

% corrosion cell
scan_path = [raw_path 'mah_08_corrosion_cell_B'];
rot_axis_offset = -2 / raw_bin;
rot_axis_tilt = -0.003;
ADD

% TODO
scan_path = [raw_path 'mah_10_13R_top'];
excentric_rot_axis = 1;
rot_axis_offset = 1078 / raw_bin;
rot_axis_tilt = -0.00267; % about -.15 degrees
ADD

scan_path = [ raw_path '/mah_10_13R_bottom'];
rot_axis_offset = 539;
rot_axis_tilt = -0.0023;
ADD

scan_path = [ raw_path 'mah_11_20R_top'];
rot_axis_offset = 538.5;
rot_axis_tilt = -0.0024;
ADD

scan_path = [ raw_path 'mah_11_20R_bottom'];
rot_axis_offset = 538.5;
rot_axis_tilt = -0.0024;
ADD('r')

scan_path = [ raw_path 'mah_15_57R'];
rot_axis_offset = -2.5 / raw_bin;
rot_axis_tilt = -0.0028;
ADD

scan_path = [ raw_path 'mah_15_57R'];
rot_axis_offset = -2.5 / raw_bin;
rot_axis_tilt = -0.0028;
do_phase_retrieval = 1;
ADD('r')

scan_path = [ raw_path 'mah_16_57R_load'];
rot_axis_offset = -2.5;
rot_axis_tilt = -0.0028;
ADD

scan_path = [ raw_path 'mah_17_57R_load_middle'];
rot_axis_offset = 88.25;
rot_axis_tilt = -0.0025;
ADD

scan_path = [ raw_path 'mah_18_57R_load_top'];
rot_axis_offset = 88.25;
rot_axis_tilt = -0.003;
ADD

% 3.4.17
scan_path = [ raw_path 'mah_20_4L_bottom'];
rot_axis_offset = 12.5;
rot_axis_tilt = -0.0029;
ADD

% movement artefacts
scan_path = [ raw_path 'mah_20_4L_top'];
rot_axis_offset = 12.5;
rot_axis_tilt = -0.0029;
ADD

do_phase_retrieval = 1;
ADD('r')

scan_path = [ raw_path 'mah_22_50L_top'];
rot_axis_offset = 88.5;
rot_axis_tilt = -0.0028;
ADD

scan_path = [ raw_path 'mah_23_50L_top'];
rot_axis_offset = 5;
rot_axis_tilt = -0.004160;
ADD

scan_path = [ raw_path 'mah_24_50L_top_load'];
rot_axis_offset = 4.5;
rot_axis_tilt = -0.004;
do_phase_retrieval = 1;
ADD('r')

scan_path = [ raw_path 'mah_28_15R_top'];
rot_axis_offset = 8.25;
rot_axis_tilt = -0.003;
ADD

scan_path = [ raw_path 'mah_29_15R_top_occd125_withpaper'];
rot_axis_offset = -3.5 / raw_bin;
rot_axis_tilt = -0.003;
ADD

scan_path = [ raw_path 'mah_29_15R_top_occd125_withpaper'];
rot_axis_offset = -3.5 / raw_bin;
rot_axis_tilt = -0.003;
do_phase_retrieval = 1;
ADD('r')

scan_path = [ raw_path 'mah_30_15R_top_occd125_withoutpaper'];
rot_axis_offset = -3.5 / raw_bin;
rot_axis_tilt = -0.003;
ADD

scan_path = [ raw_path 'mah_30_15R_top_occd125_withoutpaper'];
rot_axis_offset = -3.5 / raw_bin;
rot_axis_tilt = -0.003;
do_phase_retrieval = 1;
ADD('r')

scan_path = [ raw_path 'mah_32_15R_top_occd800_withoutpaper'];
rot_axis_offset = -79.0 / raw_bin;
rot_axis_tilt = -0.002;
do_phase_retrieval = 1;
ADD

scan_path = [ raw_path 'mah_33_50L_occd400_bottom'];
rot_axis_offset = -40 / raw_bin;
rot_axis_tilt = 0.00137;
do_phase_retrieval = 1;
ADD

scan_path = [ raw_path 'mah_33_50L_occd400_top'];
rot_axis_offset = 5 / raw_bin;
rot_axis_tilt = 0.00135;
do_phase_retrieval = 1;
ADD

scan_path = [ raw_path 'mah_35_1R_bottom'];
rot_axis_offset = 2 / raw_bin;
rot_axis_tilt = 0.00158;
do_phase_retrieval = 0;
ADD

scan_path = [ raw_path 'mah_36_1R_top'];
rot_axis_offset = 2 / raw_bin;
rot_axis_tilt = 0.00158;
ADD

scan_path = [ raw_path 'mah_37_10R_bottom'];
raw_roi = [1 1100];
rot_axis_offset = 0.6 / raw_bin;
rot_axis_tilt = 0.0015;
ADD

scan_path = [ raw_path 'mah_38_10R_top'];
rot_axis_offset = 0.6 / raw_bin;
rot_axis_tilt = 0.0015;
ADD

scan_path = [ raw_path 'mah_39_3L_bottom'];
ADD

scan_path = [ raw_path 'mah_40_3L_top'];
ADD

scan_path = [ raw_path 'mah_41_9R_bottom'];
ADD

scan_path = [ raw_path 'mah_42_9R_top'];
ADD

% Straw: no proper reco possible due to movment

% corroded screw: movement
scan_path = [raw_path 'mah_straw_01'];
rot_axis_offset = 2 / raw_bin;
rot_axis_tilt = -0.003;
ADD

% corroded screw: movement
scan_path = [raw_path 'mah_straw_02'];
rot_axis_offset = -1.5 / raw_bin;
rot_axis_tilt = -0.003;
ADD

% corroded screw
scan_path = [raw_path 'mah_straw_03'];
rot_axis_offset = -1.25 / raw_bin;
rot_axis_tilt = -0.0027;
% time-varying bright spots: for nn=1:40,imsc(flat(0+(1:400),0+(1:200),nn)',[000 9000]),pause(1),end

% corroded screw
scan_path = [raw_path 'mah_straw_04'];
rot_axis_offset = -1.75 / raw_bin;
rot_axis_tilt = -0.003;
ADD

% corroded screw
scan_path = [raw_path 'mah_straw_05'];
rot_axis_offset = -2.5 / raw_bin;
rot_axis_tilt = -0.0025;

% corroded screw
scan_path = [raw_path 'mah_straw_06'];
rot_axis_offset = -0 / raw_bin;
rot_axis_tilt = -0.003;

% corroded screw
scan_path = [raw_path 'mah_straw_2_00'];
rot_axis_offset = -0.4 / raw_bin;
rot_axis_tilt = -0.003;
write_sino = 1; 

% corroded screw
scan_path = [raw_path 'mah_straw_2_01'];
rot_axis_offset = -0 / raw_bin;
rot_axis_tilt = -0.003;
write_sino = 1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p05_reco_loop( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
