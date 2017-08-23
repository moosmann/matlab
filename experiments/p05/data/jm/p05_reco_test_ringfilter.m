function p05_reco_test_ringfilter( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
if nargin < 1
    SUBSETS = [];
end
if nargin < 2
    RUN_RECO = 0;
end
if nargin < 3
    PRINT_PARAMETERS = 'parfolder';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFAULT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scan_path = '';
raw_roi = [299 1492];
scan_path = '';
raw_bin = 1;
excentric_rot_axis = 0;
crop_at_rot_axis = 0;
stitch_projections = 0; 
proj_range = 1; 
ref_range = 1; 
correlation_method =  'ssim-ml';
do_phase_retrieval = 0;
do_tomo = 1;
fbp_filter_padding = 1;
fbp_filter_type = 'Ram-Lak';'linear';
fpb_filter_freq_cutoff = 1; 
ring_filter = 1;
ring_filter_method = 'wavelet-fft';
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
write_8bit_binned = 0;
write_16bit = 0; 
subfolder_reco = '';
gpu_ind = 1;

ADD('d')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETER / DATA SETS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scan_path = '/asap3/petra3/gpfs/p05/2017/data/11002845/raw/ste_03_g4l_aa';
parfolder = 'jm';
out_path = '/asap3/petra3/gpfs/p05/2017/data/11002845/processed/ste_03_g4l';
rot_axis_offset = -3.6/ raw_bin;
flat_corr_area1 = [1 floor(100/raw_bin)]; 
flat_corr_area2 = [0.1 0.7]; %correlation area: index vector or relative/absolute position of [first pix, last pix]
ADD

%% 2016 nov
raw = '/asap3/petra3/gpfs/p05/2016/data/11001978/raw/';
scan_path = [raw 'mah_01'];
raw_roi = [293 1492];
rot_axis_offset = -135.5 / raw_bin;
rot_axis_tilt = -0.003;

dec_levels = 6;
wname = 'db30';
parfolder = sprintf( 'test_ringfilt_%s_%s_decNum%u', ring_filter_method, wname, dec_levels );
ADD

ring_filter_method = 'jm';
parfolder = ['test_ringfilt_' ring_filter_method];
ADD

ring_filter = 0;        
ring_filter_method = 'none';
parfolder = ['test_ringfilt_' ring_filter_method];
ADD

%% 2017 may
raw = '/asap3/petra3/gpfs/p05/2017/data/11003950/raw/';

% Combined wavelet FFT
scan_path = [raw 'syn13_55L_Mg10Gd_12w_load_00'];
ring_filter_method = 'wavelet-fft';
wname = 'db25';
sigma = 2.4;
dec_levels = 1;
parfolder = sprintf('test_ringfilt_wavelet-fft_%s_decLev%s', wname, sprintf('%u', dec_levels) );
ADD

dec_levels = dec_levels + 1;
parfolder = sprintf('test_ringfilt_wavelet-fft_decNum%u', dec_levels);
ADD

dec_levels = dec_levels + 1;
parfolder = sprintf('test_ringfilt_wavelet-fft_decNum%u', dec_levels);
ADD
 
dec_levels = dec_levels + 1;
parfolder = sprintf('test_ringfilt_wavelet-fft_decNum%u', dec_levels);
ADD

dec_levels = dec_levels + 1;
parfolder = sprintf('test_ringfilt_wavelet-fft_decNum%u', dec_levels);
ADD

dec_levels = dec_levels + 1;
parfolder = sprintf('test_ringfilt_wavelet-fft_decNum%u', dec_levels);
ADD

dec_levels = dec_levels + 1;
parfolder = sprintf('test_ringfilt_wavelet-fft_decNum%u', dec_levels);
ADD


% None
ring_filter = 0;
parfolder = 'test_ringfilt_none';
ADD

% Simple JM
ring_filter = 1;
ring_filter_method = 'jm';
parfolder = 'test_ringfilt_jm';
ADD

% DB30
ring_filter = 1;
ring_filter_method = 'wavelet-fft';
wname = 'db30'; sigma = 2.4; dec_levels = 3;
parfolder = sprintf('test_ringfilt_wavelet-fft_%s_decNum%u', wname, dec_levels);
ADD

dec_levels = dec_levels + 1;
parfolder = sprintf('test_ringfilt_wavelet-fft_%s_decNum%u', wname, dec_levels);
ADD

dec_levels = dec_levels + 1;
parfolder = sprintf('test_ringfilt_wavelet-fft_%s_decNum%u', wname, dec_levels);
ADD

dec_levels = dec_levels + 1;
parfolder = sprintf('test_ringfilt_wavelet-fft_%s_decNum%u', wname, dec_levels);
ADD

dec_levels = dec_levels + 1;
parfolder = sprintf('test_ringfilt_wavelet-fft_%s_decNum%u', wname, dec_levels);
ADD(1)

%% Synload
scan_path = '/asap3/petra3/gpfs/p05/2016/commissioning/c20160913_000_synload/raw/mg5gd_02_1w';
ring_filter = 1;
%raw_roi = [141 1940];
raw_roi = [701 1100];
ring_filter_method = 'jm';
parfolder = 'test_ringfilt_jm';
ADD

ring_filter_method = 'wavelet-fft';
wname = 'db15'; sigma = 2.4; dec_levels = 1:3;
parfolder = sprintf('test_ringfilt_wavelet-fft_%s_decLev%s', wname, sprintf('%u', dec_levels) );
ADD

wname = 'db15'; sigma = 2.4; dec_levels = 1:7;
parfolder = sprintf('test_ringfilt_wavelet-fft_%s_decLev%s', wname, sprintf('%u', dec_levels) );
ADD

wname = 'db25'; sigma = 2.4; dec_levels = 1:3;
parfolder = sprintf('test_ringfilt_wavelet-fft_%s_decLev%s', wname, sprintf('%u', dec_levels) );
ADD

wname = 'db25'; sigma = 2.4; dec_levels = 1:7;
parfolder = sprintf('test_ringfilt_wavelet-fft_%s_decLev%s', wname, sprintf('%u', dec_levels) );
ADD

wname = 'db30'; sigma = 2.4; dec_levels = 1:3;
parfolder = sprintf('test_ringfilt_wavelet-fft_%s_decLev%s', wname, sprintf('%u', dec_levels) );
ADD

wname = 'db30'; sigma = 2.4; dec_levels = 1:7;
parfolder = sprintf('test_ringfilt_wavelet-fft_%s_decLev%s', wname, sprintf('%u', dec_levels) );
ADD

wname = 'bior2.8'; sigma = 2.4; dec_levels = 1:3;
parfolder = sprintf('test_ringfilt_wavelet-fft_%s_decLev%s', wname, sprintf('%u', dec_levels) );
ADD

wname = 'bior2.8'; sigma = 2.4; dec_levels = 1:7;
parfolder = sprintf('test_ringfilt_wavelet-fft_%s_decLev%s', wname, sprintf('%u', dec_levels) );
ADD

wname = 'db30'; sigma = 2.4; dec_levels = 1:7;
fbp_filter_type = 'linear'; fpb_filter_freq_cutoff = 0.5;
parfolder = sprintf('test_ringfilt_wavelet-fft_%s_decLev%s_fbpFilt%s_fbpcut%02u', wname, sprintf('%u', dec_levels), fbp_filter_type, fpb_filter_freq_cutoff*10 );
ADD

wname = 'db30'; sigma = 2.4; dec_levels = 1:7;
fbp_filter_type = 'linear'; fpb_filter_freq_cutoff = 1;
parfolder = sprintf('test_ringfilt_wavelet-fft_%s_decLev%s_fbpFilt%s_fbpcut%02u', wname, sprintf('%u', dec_levels), fbp_filter_type, fpb_filter_freq_cutoff*10 );
ADD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p05_reco_loop( SUBSETS, RUN_RECO, PRINT_PARAMETERS)
