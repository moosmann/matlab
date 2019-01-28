function p05_reco_hnee( SUBSETS, RUN_RECO, PRINT_PARAMETERS)

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
raw_bin = 2;
excentric_rot_axis = 0;
crop_at_rot_axis = 0;
stitch_projections = 0; 
proj_range = 1; 
ref_range = 1; 
correlation_method =  'ssim-ml';
do_phase_retrieval = 0;
phase_retrieval_method = 'tie';'qp';'qpcut';
phase_retrieval_reg_par = 2.5; 
phase_retrieval_bin_filt = 0.15; 
phase_retrieval_cutoff_frequ = 1 * pi; 
phase_padding = 1; 
do_tomo = 1;
butterworth_filter = 0; 
butterworth_order = 1;
butterworth_cutoff_frequ = 0.5;
fbp_filter_padding = 1;
fbp_filter_freq_cutoff = 1;
ring_filter = 1;
dec_levels = 7;
wname = 'db25';
sigma = 2.4;
rot_axis_offset = [];
rot_axis_tilt = 0;
parfolder = '';
write_to_scratch = 1;
write_flatcor = 0;
write_phase_map = 0; 
write_sino = 0; 
write_sino_phase = 0; 
write_reco = 1; 
write_float = 1; 
write_float_binned = 1; 
write_8bit = 0;
write_8bit_binned = 0;
write_16bit = 0; 

gpu_ind = [];

% Set default. Allows parameters to be changed before first data set is added
ADD('default') % or use SET_DEFAULT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETER / DATA SETS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

raw_path = '/asap3/petra3/gpfs/p05/2017/data/11003063/raw/';

scan_path = [raw_path 'hnee_01_hw_hk776_bn161514'];
ADD

scan_path = [raw_path 'hnee_02_hw_hk776_bn161514_d300'];
ADD

%% No butterworth filter
raw_roi = [101 2100];
raw_bin = 1;
phase_bin = 1; 
butterworth_filter = 0; 
butterworth_order = 1;
butterworth_cutoff_frequ = 0.5;
subfolder_reco = 'noBWfilt';

scan_path = [raw_path 'hnee_03_hw_hk776_bn161514_d150'];
subfolder_reco = '';
rot_axis_offset = 9 / raw_bin;
do_phase_retrieval = 1;
phase_retrieval_method = 'tie';
phase_retrieval_reg_par = 2.5; 
ADD

phase_retrieval_method = 'qp';
phase_retrieval_bin_filt = 0.1; 
ADD

phase_retrieval_method = 'qpcut';
phase_retrieval_cutoff_frequ = 1 * pi;
ADD

phase_retrieval_method = 'qpcut';
phase_retrieval_cutoff_frequ = 2 * pi;
ADD

raw_roi = [201 2100];
raw_bin = 2;
scan_path = [raw_path 'hnee_04_hw_hk776_bn1513_1016b'];
do_phase_retrieval = 1;
phase_retrieval_method = 'tie';
phase_retrieval_reg_par = 2.5; 
rot_axis_offset = 2*8.75 / raw_bin;
correlation_method =  'ssim-ml';
parfolder = correlation_method;
ADD

correlation_method =  'ssim';
parfolder = correlation_method;
ADD

correlation_method =  'ssim-ml';
ring_filter_method = 'wavelet-fft'; 
dec_levels = 1:3;
wname = 'db25';
sigma = 2.4;
parfolder = sprintf('ringfilt_wavelet-fft_%s_decLev%s', wname, sprintf('%u', dec_levels) );
ADD

dec_levels = 1:5;
parfolder = sprintf('ringfilt_wavelet-fft_%s_decLev%s', wname, sprintf('%u', dec_levels) );
ADD

dec_levels = 1:7;
parfolder = sprintf('ringfilt_wavelet-fft_%s_decLev%s', wname, sprintf('%u', dec_levels) );
ADD

dec_levels = 2:6;
parfolder = sprintf('ringfilt_wavelet-fft_%s_decLev%s', wname, sprintf('%u', dec_levels) );
ADD

parfolder = '';
do_phase_retrieval = 0;
ADD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p05_reco_loop( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
