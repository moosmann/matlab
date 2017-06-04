function p05_reco_test_fbp( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
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

raw_roi = [];
scan_path = '';
raw_bin = 1;
excentric_rot_axis = 0;
crop_at_rot_axis = 0;
stitch_projections = 0; 
proj_range = 1; 
ref_range = 1; 
correlation_method =  'ssim-ml';
ring_filter = 1; % ring artifact filter
ring_filter_method = 'jm';'wavelet-fft'; 
ring_filter_median_width = 11;
do_phase_retrieval = 0;
phase_retrieval_method = 'tie';'qp';'qpcut';
phase_retrieval_reg_par = 2.5; 
phase_retrieval_bin_filt = 0.15; 
phase_retrieval_cutoff_frequ = 1 * pi; 
phase_padding = 1; 
fbp_filter_type = 'Ram-Lak';'linear';
fpb_filter_freq_cutoff = 1;
fbp_filter_padding = 1;
fbp_filter_padding_method = 'symmetric';
do_tomo = 1;
ring_filter = 1;
rot_axis_offset = [];
rot_axis_tilt = [];
butterworth_filter = 1; % use butterworth filter in addition to FBP filter
butterworth_order = 1;
butterworth_cutoff_frequ = 0.5;
astra_pixel_size = 1; % size of a detector pixel: if different from one 'vol_size' needs to be ajusted
take_neg_log = [];
parfolder = '';
write_to_scratch = 0;
write_flatcor = 0;
write_phase_map = 0; 
write_sino = 0; 
write_sino_phase = 0; 
write_reco = 1; 
write_float = 0; 
write_float_binned = 1; 
write_8bit = 0;
write_8bit_binned = 0;
write_16bit = 0; 

ADD('d')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETER / DATA SETS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%sprintf('fbpFilt%s_ringFilt%uMedWid%u_bwFilt%ubwCutoff%u_phasePad%u_freqCutoff%2.0f_fbpPad%u', fbp_filter_type, ring_filter, ring_filter_median_width, butterworth_filter, 100*butterworth_cutoff_frequ, phase_padding, fpb_freq_cutoff*100, fbp_filter_padding); % parent folder to 'reco'

sprintf('fbpFilt%s_ringFilt%uMedWid%u_bwFilt%ubwCutoff%u_phasePad%u_freqCutoff%2.0f_fbpPad%u', fbp_filter_type, ring_filter, ring_filter_median_width, butterworth_filter, 100*butterworth_cutoff_frequ, phase_padding, fpb_filter_freq_cutoff*100, fbp_filter_padding); 
out_path = '/gpfs/petra3/scratch/moosmanj';
out_path = '/asap3/petra3/gpfs/p05/2016/data/11001978/scratch_cc/c20160803_001_pc_test';

raw_path = '/asap3/petra3/gpfs/p05/2017/data/11003950/raw/';
raw_roi = [265 1464];

scan_path = [raw_path 'syn01_48L_PEEK_12w_b'];
do_phase_retrieval = 1;phase_retrieval_method = 'tie';phase_retrieval_reg_par = 2.5; 
ADD

scan_path = [raw_path 'syn01_48L_PEEK_12w_b'];
parfolder = 'fbp_fft_wo_symmetric_option';
ADD

scan_path = [raw_path 'syn01_48L_PEEK_12w_b'];
do_phase_retrieval = 1;phase_retrieval_method = 'tie';phase_retrieval_reg_par = 2.5; 
parfolder = 'fbp_fft_wo_symmetric_option';
ADD

scan_path = [raw_path 'syn01_48L_PEEK_12w_b'];
do_phase_retrieval = 1;phase_retrieval_method = 'tie';phase_retrieval_reg_par = 2.5; 
parfolder = 'fbp_fft_wo_symmetric_option__nopadding';
phase_padding = 0;
fbp_filter_padding = 0;
ADD

scan_path = [raw_path 'syn01_48L_PEEK_12w_b'];
parfolder = 'fbp_fft_with_symmetric_option';
ADD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p05_reco_loop( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
