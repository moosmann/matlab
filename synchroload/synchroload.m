% Synchroload project related scripts

%% Reco loop scripts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
edit synchroload2016nov_11001978.m
edit synchroload2016sept_c20160913.m
edit synchroload2016sept_c20160913_check.m.m % c20160920_000_diana c20160913_000_synload
edit synchroload2017may_11003950.m
edit synchroload2017nov03_11003773
edit synchroload2017nov13_11004016.m
edit synchroload2017nov23_11003288.m
edit synchroload2017oct_11003440.m
edit synchroload2018apr_11004263.m
edit synchroload2018may_11004936.m
edit synchroload2018nov_11005553.m
edit synchroload2019may_11005842.m
edit synchroload2019july_11006704.m
edit synchroload2019dec_11006991.m

%% Other scripts
edit synchroload_scans.m
edit synchroload_figureMeetingApr2017.m
edit synchroload_radiography.m
edit synchroload_renorm_slices.m
edit synchroload_resample_data_for_segmentation_old.m
edit synchroload_resample_data_for_segmentation_autodetect.m
edit synchroload_resample_data_for_segmentation_listmode.m
edit synchroload_sequ_force.m
edit synchroload_sequ.m
edit synchroload_syn13.m
edit synchroload_syn166.m
edit stitch_volumes.m
edit dose_synchroload.m
edit dose_bone_test.m

%% Segmentation
% High quality
'/asap3/petra3/gpfs/external/2019/data/50000258/processed/resampled/113734/segmentation_hanna';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load force value: create figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.steps = [];
p.out_path = '';'/gpfs/petra3/scratch/moosmanj/nextcloud';

%% 2019 July 11006704 LTP IV
p.raw_path = '/asap3/petra3/gpfs/p05/2019/data/11006704/raw';

p.steps = [];
p.out_path = '';
p.readhdf5 = 1;
% load_cell_20N_large
%calibration factor V-to-N:       2.45250 N/V
p.adc2force = 2.45250;

%200N load cell large
% Load sequence 2019-07-13 18:00
%load_sequence, 'syn004_71L_Mg5Gd_12w'
%calibration factor V-to-N:       8.35605 N/V
p.adc2force = 8.35605;

%p.scan_name = 'syn001_35R_Ti_8w';
%p.scan_name = 'syn002_65L_Mg5Gd_12w';
p.steps = [6:21 26:33];
p.scan_name = 'syn004_71L_Mg5Gd_12w';
%p05_load_force_values( p )
p05_load_force_values( p )

%load_cell_200N_large, weight, voltage
weight = 5.31;
voltage = 2.65;
gconst = 9.81;
p.adc2force = weight * gconst / voltage;
p.scan_name = 'syn007_47L_Peek_12w';
p.steps = [5:22 24 26:33];
p.scan_name = 'syn008_47R_Ti_12w';


%% 2016 

%% 2017 May 11003950
p.raw_path = '/asap3/petra3/gpfs/p05/2017/data/11003950/raw';
p.scan_name = 'syn13_55L_Mg10Gd_12w_load';
p.steps = [];
p.out_path = '';
p.adc2force = 19.9966; % Not needed for plotting
p.readhdf5 = 0;
p05_load_force_values( p )

%% 2017 Nov 13 11004016
p.raw_path = '/asap3/petra3/gpfs/p05/2017/data/11004016/raw';
p.steps = [];
p.out_path = '';
p.adc2force = 5.31 / 10.62*9.81; % Not needed for plotting
p.readhdf5 = 1;
p.scan_name = 'syn002_6L_PEEK_4w';
p.scan_name = 'syn003_92L_Mg10Gd_4w';
p.scan_name = 'syn004_84L_Mg10Gd_4w';
p.scan_name = 'syn005_81L_Mg5Gd_8w';
p.scan_name = 'syn006_75R_Mg10Gd_8w';
p.scan_name = 'syn007_94L_Mg10Gd_8w';
p.scan_name = 'syn008_76R_Mg10Gd_8w';
p.scan_name = 'syn009_32R_PEEK_8w';
p.scan_name = 'syn010_19R_PEEK_4w';
p.scan_name = 'syn011_14R_PEEK_4w';
p.scan_name = 'syn012_79L_Mg5Gd_8w';
%p05_load_force_values( p )
s = {
'syn002_6L_PEEK_4w';
'syn003_92L_Mg10Gd_4w';
'syn004_84L_Mg10Gd_4w';
'syn005_81L_Mg5Gd_8w';
'syn006_75R_Mg10Gd_8w';
'syn007_94L_Mg10Gd_8w';
'syn008_76R_Mg10Gd_8w';
'syn009_32R_PEEK_8w';
'syn010_19R_PEEK_4w';
'syn011_14R_PEEK_4w';
'syn012_79L_Mg5Gd_8w';};
for scan_name = s'
    p.scan_name = scan_name{1};
    p05_load_force_values( p )
end

%% 2017 Nov 23 11003288
p.raw_path = '/asap3/petra3/gpfs/p05/2017/data/11003288/raw';
p.steps = [];
p.out_path = '';
p.adc2force = 4.81879;
p.scan_name = 'syn166_104R_Mg10Gd_4w';
p.readhdf5 = 0;
p05_load_force_values( p )

%% 2018 April/May 11004263
% None

%% 2018 May 25 11004936
p.raw_path = '/asap3/petra3/gpfs/p05/2018/data/11004936/raw';

p.scan_name = 'syn007_56R_Mg5Gd_12w';
p05_load_force_values( p )

p.scan_name = 'syn005_55R_Mg5Gd_12w_load';
p05_load_force_values( p.raw_path, p.scan_name, p.steps, p.out_path )

%% 2018 Nov 11005553 LTP III 
p.raw_path = '/asap3/petra3/gpfs/p05/2018/data/11005553/raw';

p.scan_name = 'syn004_24R_PEEK_8w';
p05_load_force_values( p.raw_path, p.scan_name)

p.scan_name = 'syn008_24R_PEEK_8w';
p05_load_force_values( p.raw_path, p.scan_name)
 
p.scan_name = 'syn013_105L_Mg5Gd_4w';
p05_load_force_values( p.raw_path, p.scan_name,[],'',10 )

p.scan_name = 'syn014_105R_Mg10Gd_4w';
p05_load_force_values( p.raw_path, p.scan_name,[],'',10 )

p.scan_name = 'syn018_86L_Mg10Gd_4w_afterdrilling_push';
p05_load_force_values( p.raw_path, p.scan_name,[],'',10)

p.scan_name = 'syn026_femur_55L';
p05_load_force_values( p.raw_path, p.scan_name,[],'',10)

% stability test
filename = '/asap3/petra3/gpfs/p05/2018/data/11005553/raw/test/easyform_15x2N_relax60min_nexus.h5';
p.out_path = '';
p.adc2force = 0;
p05_load_plot_from_hdf5( filename, p.out_path, p.adc2force );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load sequence animation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2019 July 11006704
p.regdir = 'x';
p.steps = [];
p.register = 0;
p.outlier_thresh = 0.01;
p.proc_path = '/asap3/petra3/gpfs/p05/2019/data/11006704/processed';
p.voxel_size = 6.4e-6;
p.barcol = 'white';
p.reco_sub = 'reco/float_rawBin5';
p.scan_name = 'syn004_71L_Mg5Gd_12w';
%p.scan_name = 'syn007_47L_Peek_12w';
%p.scan_name = 'syn008_47R_Ti_12w';
%vol = p05_load_sequ( p );

%% 2017 Nov 13 11004016
p.proc_path = '/asap3/petra3/gpfs/p05/2017/data/11004016/processed';
%p.scan_name = 'syn002_6L_PEEK_4w';
%p.scan_name = 'syn003_92L_Mg10Gd_4w';
%p.scan_name = 'syn004_84L_Mg10Gd_4w'; p.regdir = 'y';
%p.scan_name = 'syn005_81L_Mg5Gd_8w';
%p.scan_name = 'syn006_75R_Mg10Gd_8w';
%p.scan_name = 'syn007_94L_Mg10Gd_8w';
%p.scan_name = 'syn008_76R_Mg10Gd_8w';
%p.scan_name = 'syn009_32R_PEEK_8w';
%p.scan_name = 'syn010_19R_PEEK_4w';
%p.scan_name = 'syn011_14R_PEEK_4w'; p.regdir = 'y';
%p.scan_name = 'syn012_79L_Mg5Gd_8w'; % empty scans [1:3 14:20];

%% 2018 11004936

p.proc_path = '/asap3/petra3/gpfs/p05/2018/data/11004936/processed';
p.reco_sub = 'reco/float_rawBin2_phaseBin2';
'syn007_56R_Mg5Gd_12w_';
p.scan_name = 'syn005_55R_Mg5Gd_12w_load';
p.regdir = 'x';
p.steps = [];
p.register = 1;
p.outlier_thresh = 0.01;
[vol, vol_reg] = p05_load_sequ( p );

%% 2018 11005553

p.proc_path = '/asap3/petra3/gpfs/p05/2018/data/11005553/processed';
p.reco_sub = 'reco/float_rawBin2';
p.scan_name = 'syn004_24R_PEEK_8w';
p.regdir = 'x';
p.steps = [];
p.register = 0;
p.outlier_thresh = 0.01;
[vol, vol_reg] = p05_load_sequ( p );


p.proc_path = '/asap3/petra3/gpfs/p05/2018/data/11005553/processed';
p.reco_sub = 'reco/float_rawBin2';
p.scan_name = 'syn013_105L_Mg5Gd_4w';
p.regdir = 'x';
p.steps = [];
p.register = 1;
p.outlier_thresh = 0.01;
[vol, vol_reg] = p05_load_sequ( p );

p.proc_path = '/asap3/petra3/gpfs/p05/2018/data/11005553/processed';
p.reco_sub = 'reco/float_rawBin3';
p.scan_name = 'syn018_86L_Mg10Gd_4w_afterdrilling_push';
p.regdir = 'x';
p.steps = [];
p.register = 1;
p.outlier_thresh = 0.01;
[vol, vol_reg] = p05_load_sequ( p );

%%
p.proc_path = '/asap3/petra3/gpfs/p05/2018/data/11005553/processed';
p.reco_sub = 'reco/float_rawBin3_phaseBin2';
p.scan_name = 'syn014_105R_Mg10Gd_4w';
p.regdir = 'x';
p.steps = [];
p.register = 0;
p.outlier_thresh = 0.01;
[vol, vol_reg] = p05_load_sequ( p );

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
p.out_path = '/asap3/petra3/gpfs/p05/2017/data/11003288/processed/syn166_104R_Mg10Gd_4w/dvc/HxUniformVectorField3';
out_name = 'vec';
steps_to_process = [];
write_single_component_vector_binaries_to_HxUniformVectorField3( input_folder, folder_pattern, filename_pattern, shape, dtype, machinefmt, p.out_path, out_name, steps_to_process );

%% syn13 4D volume cropped, DVC with blockiness
input_folder = '/asap3/petra3/gpfs/p05/2017/data/11003950/processed/syn13_55L_Mg10Gd_12w_load/dvc/raw';
folder_pattern = 'out*';
filename_pattern = '*vec*.raw';
shape = [700 600 500];
dtype = 'single';
machinefmt = 'b';
p.out_path = '/asap3/petra3/gpfs/p05/2017/data/11003950/processed/syn13_55L_Mg10Gd_12w_load/dvc/HxUniformVectorField3';
out_name = 'vec';
steps_to_process = [];
write_single_component_vector_binaries_to_HxUniformVectorField3( input_folder, folder_pattern, filename_pattern, shape, dtype, machinefmt, p.out_path, out_name, steps_to_process );

%% DVC Divergence

%% syn13
input_folder = '/asap3/petra3/gpfs/p05/2017/data/11003950/processed/syn13_55L_Mg10Gd_12w_load/dvc/raw';
folder_pattern = 'out*';
filename_pattern = 'stk*.raw';
shape = [700 600 500];
dtype = 'single';
machinefmt = 'b';
p.out_path = '/asap3/petra3/gpfs/p05/2017/data/11003950/processed/syn13_55L_Mg10Gd_12w_load/dvc/HxUniformScalarField';
out_name = 'div';
steps_to_process = [];
write_divergence_vector_binaries_to_HxUniformScalarField( input_folder, folder_pattern, filename_pattern, shape, dtype, machinefmt, p.out_path, out_name, steps_to_process )

%% syn166
input_folder = '/asap3/petra3/gpfs/p05/2017/data/11003288/processed/syn166_104R_Mg10Gd_4w/dvc';
folder_pattern = 'out*';
filename_pattern = 'vec*.vol';
shape = [1240 1120 20];
dtype = 'single';
machinefmt = 'l';
p.out_path = '/asap3/petra3/gpfs/p05/2017/data/11003288/processed/syn166_104R_Mg10Gd_4w/dvc/HxUniformScalarField';
out_name = 'div';
steps_to_process = [];
write_divergence_vector_binaries_to_HxUniformScalarField( input_folder, folder_pattern, filename_pattern, shape, dtype, machinefmt, p.out_path, out_name, steps_to_process )

%% Dose %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dose_synchroload

%% Segmentation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

synchroload_resample_data_for_segmentation

%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
