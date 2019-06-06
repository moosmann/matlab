% Synchroload project related scripts

%% Reco loop scripts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
edit synchroload2016nov_11001978.m
edit synchroload2016sept_check.m
edit synchroload2016sept.m % c20160920_000_diana c20160913_000_synload
edit synchroload2017may_11003950.m
edit synchroload2017nov13_11004016.m
edit synchroload2017nov23_11003288.m
edit synchroload2017nov23_11003288_oldDoubleCheck.m
edit synchroload2017oct_11003440_1.m
edit synchroload2017oct_11003440.m
edit synchroload2018apr_11004263.m
edit synchroload2018may_11004936.m
edit synchroload2018nov_11005553.m
edit synchroload2019may_11005842

%% Other scripts
edit synchroload_figureMeetingApr2017.m
edit synchroload_radiography.m
edit synchroload_renorm_slices.m
edit synchroload_resample_data_for_segmentation.m
edit synchroload_sequ_force.m
edit synchroload_sequ.m
edit synchroload_syn13.m
edit synchroload_syn166.m
edit stitch_volumes
edit find_volume_stitch_parameter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load force value: create figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
steps = [];
out_path = '';'/gpfs/petra3/scratch/moosmanj/nextcloud';

%% 2016 

%% 2017 May 11003950
raw_path = '/asap3/petra3/gpfs/p05/2017/data/11003950/raw';
scan_name = 'syn13_55L_Mg10Gd_12w_load';
steps = [];
out_path = '';
adc2force = 19.9966; % Not needed for plotting
readhdf5 = 0;
p05_load_force_values( raw_path, scan_name, steps, out_path, adc2force, readhdf5 )

%% 2017 Nov 13 11004016
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

%% 2017 Nov 23 11003288
raw_path = '/asap3/petra3/gpfs/p05/2017/data/11003288/raw';
steps = [];
out_path = '';
adc2force = 4.81879;
scan_name = 'syn166_104R_Mg10Gd_4w';
readhdf5 = 0;
p05_load_force_values( raw_path, scan_name, steps, out_path, adc2force, readhdf5 )

%% 2018 April/May 11004263
% None

%% 2018 May 25 11004936
raw_path = '/asap3/petra3/gpfs/p05/2018/data/11004936/raw';

scan_name = 'syn007_56R_Mg5Gd_12w';
p05_load_force_values( raw_path, scan_name, steps, out_path )

scan_name = 'syn005_55R_Mg5Gd_12w_load';
p05_load_force_values( raw_path, scan_name, steps, out_path )

%% 2018 Nov 11005553 LTP III 
raw_path = '/asap3/petra3/gpfs/p05/2018/data/11005553/raw';

scan_name = 'syn004_24R_PEEK_8w';
p05_load_force_values( raw_path, scan_name)

scan_name = 'syn008_24R_PEEK_8w';
p05_load_force_values( raw_path, scan_name)
 
scan_name = 'syn013_105L_Mg5Gd_4w';
p05_load_force_values( raw_path, scan_name,[],'',10 )

scan_name = 'syn014_105R_Mg10Gd_4w';
p05_load_force_values( raw_path, scan_name,[],'',10 )

scan_name = 'syn018_86L_Mg10Gd_4w_afterdrilling_push';
p05_load_force_values( raw_path, scan_name,[],'',10)

scan_name = 'syn026_femur_55L';
p05_load_force_values( raw_path, scan_name,[],'',10)

% stability test
filename = '/asap3/petra3/gpfs/p05/2018/data/11005553/raw/test/easyform_15x2N_relax60min_nexus.h5';
out_path = '';
adc2force = 0;
p05_load_plot_from_hdf5( filename, out_path, adc2force );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load sequence animation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2017 Nov 13 11004016
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

%% 2018 11004936

proc_path = '/asap3/petra3/gpfs/p05/2018/data/11004936/processed';
reco_subfolder = 'reco/float_rawBin2_phaseBin2';
'syn007_56R_Mg5Gd_12w_';
scan = 'syn005_55R_Mg5Gd_12w_load';
regdir = 'x';
steps = [];
register = 1;
outlier_thresh = 0.01;
[vol, vol_reg] = p05_load_sequ(proc_path, scan, reco_subfolder, regdir, steps, outlier_thresh, register);

%% 2018 11005553

proc_path = '/asap3/petra3/gpfs/p05/2018/data/11005553/processed';
reco_subfolder = 'reco/float_rawBin2';
scan = 'syn004_24R_PEEK_8w';
regdir = 'x';
steps = [];
register = 1;
outlier_thresh = 0.01;
[vol, vol_reg] = p05_load_sequ(proc_path, scan, reco_subfolder, regdir, steps, outlier_thresh, register);

proc_path = '/asap3/petra3/gpfs/p05/2018/data/11005553/processed';
reco_subfolder = 'reco/float_rawBin2';
scan = 'syn013_105L_Mg5Gd_4w';
regdir = 'x';
steps = [];
register = 1;
outlier_thresh = 0.01;
[vol, vol_reg] = p05_load_sequ(proc_path, scan, reco_subfolder, regdir, steps, outlier_thresh, register);


proc_path = '/asap3/petra3/gpfs/p05/2018/data/11005553/processed';
reco_subfolder = 'reco/float_rawBin3';
scan = 'syn018_86L_Mg10Gd_4w_afterdrilling_push';
regdir = 'x';
steps = [];
register = 1;
outlier_thresh = 0.01;
[vol, vol_reg] = p05_load_sequ(proc_path, scan, reco_subfolder, regdir, steps, outlier_thresh, register);

proc_path = '/asap3/petra3/gpfs/p05/2018/data/11005553/processed';
reco_subfolder = 'reco/float_rawBin3_phaseBin2';
scan = 'syn014_105R_Mg10Gd_4w';
regdir = 'x';
steps = [];
register = 1;
outlier_thresh = 0.01;
[vol, vol_reg] = p05_load_sequ(proc_path, scan, reco_subfolder, regdir, steps, outlier_thresh, register);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DVC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert vector fields to Amira / Avizo compliant vector format

%% syn166
input_folder = '/asap3/petra3/gpfs/p05/2017/data/11003288/processed/syn166_104R_Mg10Gd_4w/dvc';
folder_pattern = 'out*';
filename_pattern = 'vec*.vol';
shape = [1240 1120 20];
dtype = 'single';
machinefmt = 'l';
out_path = '/asap3/petra3/gpfs/p05/2017/data/11003288/processed/syn166_104R_Mg10Gd_4w/dvc/HxUniformVectorField3';
out_name = 'vec';
steps_to_process = [];
write_single_component_vector_binaries_to_HxUniformVectorField3( input_folder, folder_pattern, filename_pattern, shape, dtype, machinefmt, out_path, out_name, steps_to_process );

%% syn13 4D volume cropped, DVC with blockiness
input_folder = '/asap3/petra3/gpfs/p05/2017/data/11003950/processed/syn13_55L_Mg10Gd_12w_load/dvc/raw';
folder_pattern = 'out*';
filename_pattern = '*vec*.raw';
shape = [700 600 500];
dtype = 'single';
machinefmt = 'b';
out_path = '/asap3/petra3/gpfs/p05/2017/data/11003950/processed/syn13_55L_Mg10Gd_12w_load/dvc/HxUniformVectorField3';
out_name = 'vec';
steps_to_process = [];
write_single_component_vector_binaries_to_HxUniformVectorField3( input_folder, folder_pattern, filename_pattern, shape, dtype, machinefmt, out_path, out_name, steps_to_process );

%% DVC Divergence

%% syn13
input_folder = '/asap3/petra3/gpfs/p05/2017/data/11003950/processed/syn13_55L_Mg10Gd_12w_load/dvc/raw';
folder_pattern = 'out*';
filename_pattern = 'stk*.raw';
shape = [700 600 500];
dtype = 'single';
machinefmt = 'b';
out_path = '/asap3/petra3/gpfs/p05/2017/data/11003950/processed/syn13_55L_Mg10Gd_12w_load/dvc/HxUniformScalarField';
out_name = 'div';
steps_to_process = [];
write_divergence_vector_binaries_to_HxUniformScalarField( input_folder, folder_pattern, filename_pattern, shape, dtype, machinefmt, out_path, out_name, steps_to_process )

%% syn166
input_folder = '/asap3/petra3/gpfs/p05/2017/data/11003288/processed/syn166_104R_Mg10Gd_4w/dvc';
folder_pattern = 'out*';
filename_pattern = 'vec*.vol';
shape = [1240 1120 20];
dtype = 'single';
machinefmt = 'l';
out_path = '/asap3/petra3/gpfs/p05/2017/data/11003288/processed/syn166_104R_Mg10Gd_4w/dvc/HxUniformScalarField';
out_name = 'div';
steps_to_process = [];
write_divergence_vector_binaries_to_HxUniformScalarField( input_folder, folder_pattern, filename_pattern, shape, dtype, machinefmt, out_path, out_name, steps_to_process )

%% Dose %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dose_synchroload

%% Segmentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

synchroload_resample_data_for_segmentation

%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
