function p05_reco_test_fbp(nums, doreco, print_field)

if nargin < 1
    nums = [];
end
if nargin < 2
    doreco = 0;
end
if nargin < 3
    print_field = '';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter sets to loop over %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sprintf('fbpFilt%s_ringFilt%uMedWid%u_bwFilt%ubwCutoff%u_phasePad%u_freqCutoff%2.0f_fbpPad%u', fbp_filter_type, ring_filter, ring_filter_median_width, butterworth_filter, 100*butterworth_cutoff_frequ, phase_padding, fpb_freq_cutoff*100, fbp_filter_padding); % parent folder to 'reco'

%% Default %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nn = 0;
default.raw_roi = [];
default.scan_path = '';
default.raw_bin = 1;
default.excentric_rot_axis = 0;
default.crop_at_rot_axis = 0;
default.stitch_projections = 0; 
default.proj_range = 1; 
default.ref_range = 1; 
default.correlation_method =  'ssim-ml';
default.do_phase_retrieval = 0;
default.phase_retrieval_method = 'tie';'qp';'qpcut';
default.phase_retrieval_reg_par = 2.5; 
default.phase_retrieval_bin_filt = 0.15; 
default.phase_retrieval_cutoff_frequ = 1 * pi; 
default.phase_padding = 1; 
default.do_tomo = 1;
default.ring_filter = 1;
default.rot_axis_offset = [];
default.rot_axis_tilt = [];
default.parfolder = '';
default.write_to_scratch = 0;
default.write_flatcor = 1;
default.write_phase_map = 0; 
default.write_sino = 0; 
default.write_sino_phase = 0; 
default.write_reco = 1; 
default.write_float = 1; 
default.write_float_binned = 1; 
default.write_8bit = 1;
default.write_8bit_binned = 1;
default.write_16bit = 0; 

%% Data sets %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nn = nn +1; para(nn) = default;
sprintf('fbpFilt%s_ringFilt%uMedWid%u_bwFilt%ubwCutoff%u_phasePad%u_freqCutoff%2.0f_fbpPad%u', fbp_filter_type, ring_filter, ring_filter_median_width, butterworth_filter, 100*butterworth_cutoff_frequ, phase_padding, fpb_freq_cutoff*100, fbp_filter_padding); 
out_path = '/gpfs/petra3/scratch/moosmanj';
out_path = '/asap3/petra3/gpfs/p05/2016/data/11001978/scratch_cc/c20160803_001_pc_test';

%% TEST SECTION
default.write_flatcor = 0;
default.write_8bit = 0;
default.write_8bit_binned = 0;
default.raw_roi = [265 1464];

nn = nn + 1;para(nn) = default;para(nn).scan_path = [raw 'syn01_48L_PEEK_12w_b'];
para(nn).do_phase_retrieval = 1;para(nn).phase_retrieval_method = 'tie';para(nn).phase_retrieval_reg_par = 2.5; 

nn = nn + 1;para(nn) = default;para(nn).scan_path = [raw 'syn01_48L_PEEK_12w_b'];
para(nn).parfolder = 'fbp_fft_wo_symmetric_option';

nn = nn + 1;para(nn) = default;para(nn).scan_path = [raw 'syn01_48L_PEEK_12w_b'];
para(nn).do_phase_retrieval = 1;para(nn).phase_retrieval_method = 'tie';para(nn).phase_retrieval_reg_par = 2.5; 
para(nn).parfolder = 'fbp_fft_wo_symmetric_option';

nn = nn + 1;para(nn) = default;para(nn).scan_path = [raw 'syn01_48L_PEEK_12w_b'];
para(nn).do_phase_retrieval = 1;para(nn).phase_retrieval_method = 'tie';para(nn).phase_retrieval_reg_par = 2.5; 
para(nn).parfolder = 'fbp_fft_wo_symmetric_option__nopadding';
para(nn).phase_padding = 0;
para(nn).fbp_filter_padding = 0;

nn = nn + 1;para(nn) = default;para(nn).scan_path = [raw 'syn01_48L_PEEK_12w_b'];
para(nn).parfolder = 'fbp_fft_with_symmetric_option';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p05_reco_loop( nums, doreco, print_field, para)
