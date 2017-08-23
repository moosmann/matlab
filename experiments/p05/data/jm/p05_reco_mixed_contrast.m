function p05_reco_mixed_contrast( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
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

visualOutput = 0;
interactive_determination_of_rot_axis = 0;
interactive_determination_of_rot_axis_tilt = 0;
raw_roi = [];
scan_path = '';
raw_bin = 1;
excentric_rot_axis = 0;
crop_at_rot_axis = 0;
stitch_projections = 0; 
proj_range = 1; 
ref_range = 1; 
correlation_method =  'ssim-ml';
do_phase_retrieval = 1;
phase_retrieval_method = 'tie';'qp';'qpcut';
phase_retrieval_reg_par = 2.5; 
phase_retrieval_bin_filt = 0.15; 
phase_retrieval_cutoff_frequ = 1 * pi; 
phase_padding = 1; 
do_tomo = 1;
fbp_filter_padding = 1;
fbp_filter_freq_cutoff = 1;
ring_filter = 1;
dec_levels = 7;
wname = 'db25';
sigma = 2.4;
rot_axis_offset = [];
rot_axis_tilt = 0.001;
parfolder = '';
write_to_scratch = 0;
write_flatcor = 0;
write_phase_map = 0; 
write_sino = 0; 
write_sino_phase = 0; 
write_reco = 1; 
write_float = 1; 
write_float_binned = 1; 
write_8bit = 0;
write_8bit_binned = 1;
write_16bit = 0; 
subfolder_reco = '';
gpu_ind = 1;

% Set default. Allows parameters to be changed before first data set is added
ADD('default') % or use SET_DEFAULT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETER / DATA SETS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scan_path = '/asap3/petra3/gpfs/p05/2016/data/11001978/raw/mah_32_15R_top_occd800_withoutpaper'; % too much fringes, not enough coherence probably, using standard phase retrieval reco looks blurry
ADD

scan_path = '/asap3/petra3/gpfs/p05/2016/commissioning/c20160803_001_pc_test/raw/phase_600';
ADD

scan_path = '/asap3/petra3/gpfs/p05/2016/commissioning/c20160803_001_pc_test/raw/phase_1000'; %rot_axis_offset=7.5
ADD

scan_path = '/asap3/petra3/gpfs/p05/2016/commissioning/c20160803_001_pc_test/raw/phase_1400'; %rot_axis_offset=19.5
ADD

scan_path = '/asap3/petra3/gpfs/p05/2016/commissioning/c20160803_001_pc_test/raw/we43_phase_030';
ADD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p05_reco_loop( SUBSETS, RUN_RECO, PRINT_PARAMETERS)


