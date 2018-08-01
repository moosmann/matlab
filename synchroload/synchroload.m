% Synchroload project related files and scripts

%% Reco loop scripts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%/asap3/petra3/gpfs/common/p05/jm/matlab/experiments/p05/data/moosmanj/p05_reco_loop_synchroload2018june.m
edit p05_reco_loop_synchroload2018june
%/asap3/petra3/gpfs/common/p05/jm/matlab/experiments/p05/data/moosmanj/p05_reco_loop_synchroload2018may_000.m
edit p05_reco_loop_synchroload2018may_000
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load force value: create figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
steps = [];
out_path = '';'/gpfs/petra3/scratch/moosmanj/nextcloud';

raw_path = '/asap3/petra3/gpfs/p05/2017/data/11004016/raw';
%scan_name = 'syn002_6L_PEEK_4w';
%scan_name = 'syn003_92L_Mg10Gd_4w';
%scan_name = 'syn004_84L_Mg10Gd_4w';
%scan_name = 'syn005_81L_Mg5Gd_8w';
%scan_name = 'syn006_75R_Mg10Gd_8w';
%scan_name = 'syn007_94L_Mg10Gd_8w';
%scan_name = 'syn008_76R_Mg10Gd_8w';
%scan_name = 'syn009_32R_PEEK_8w';
%scan_name = 'syn010_19R_PEEK_4w';
scan_name = 'syn011_14R_PEEK_4w';
%scan_name = 'syn012_79L_Mg5Gd_8w';

raw_path = '/asap3/petra3/gpfs/p05/2018/data/11004936/raw';

scan_name = 'syn007_56R_Mg5Gd_12w';
p05_load_force_values( raw_path, scan_name, steps, out_path )

scan_name = 'syn005_55R_Mg5Gd_12w_load';
p05_load_force_values( raw_path, scan_name, steps, out_path )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load sequence processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
proc_path = '/asap3/petra3/gpfs/p05/2017/data/11004016/processed';
%scan_name = 'syn002_6L_PEEK_4w';
%scan_name = 'syn003_92L_Mg10Gd_4w';
%scan_name = 'syn004_84L_Mg10Gd_4w'; regdir = 'y';
%scan_name = 'syn005_81L_Mg5Gd_8w';
%scan_name = 'syn006_75R_Mg10Gd_8w';
%scan_name = 'syn007_94L_Mg10Gd_8w';
%scan_name = 'syn008_76R_Mg10Gd_8w';
%scan_name = 'syn009_32R_PEEK_8w';
%scan_name = 'syn010_19R_PEEK_4w';
%scan_name = 'syn011_14R_PEEK_4w'; regdir = 'y';
%scan_name = 'syn012_79L_Mg5Gd_8w'; % empty scans [1:3 14:20];


proc_path = '/asap3/petra3/gpfs/p05/2018/data/11004936/processed';
reco_subfolder = 'reco/float_rawBin2_phaseBin2';
'syn007_56R_Mg5Gd_12w_';
scan = 'syn005_55R_Mg5Gd_12w_load';
regdir = 'x';
steps = [];
register = 1;
outlier_thresh = 0.01;
[vol, vol_reg] = p05_load_sequ(proc_path, scan, reco_subfolder, regdir, steps, outlier_thresh, register);


