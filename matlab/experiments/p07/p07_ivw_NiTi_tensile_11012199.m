%% Load sequence animation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.regdir = 'z';
p.steps = [];
p.register = 1;
p.outlier_thresh = 0.0001;
p.proc_path = '/asap3/petra3/gpfs/p07/2021/data/11012199/processed';
p.voxel_size = 2.55e-6;
p.barcol = 'white';
p.reco_sub = 'reco/float_rawBin2';
p.scan_name = 'ivw0015_referenzblau_1';
p.auto_roi_center = 0;
p.crop_roi = 0;
p05_load_sequ( p );

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.regdir = 'x';
p.steps = [];
p.register = 0;
p.outlier_thresh = 0.0001;
p.proc_path = '/asap3/petra3/gpfs/p07/2021/data/11012199/processed';
p.voxel_size = 2.55e-6;
p.barcol = 'white';
p.reco_sub = 'reco/float_rawBin2';
p.scan_name = 'ivw0015_referenzblau_2';
vol = p05_load_sequ( p );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.regdir = 'z';
p.steps = [];
p.register = 1;
p.outlier_thresh = 0.0001;
p.proc_path = '/asap3/petra3/gpfs/p07/2021/data/11012199/processed';
p.voxel_size = 2.55e-6;
p.barcol = 'white';
p.reco_sub = 'reco/float_rawBin2';
p.scan_name = 'ivw0015_referenzblau_2';
p05_load_sequ( p );

%p.scan_name = 'ivw0017_Struktur1_gruen_1'; p05_load_sequ( p );
p.scan_name = 'ivw0020_Struktur1_gruen_2'; p05_load_sequ( p );
p.scan_name = 'ivw0021_Struktur1_gruen_2b';p05_load_sequ( p );
p.scan_name = 'ivw0025_Struktur1_gruen_3';p05_load_sequ( p );
p.scan_name = 'ivw0027_Struktur2_pink_1b'; p05_load_sequ( p );
p.scan_name = 'ivw0027_Struktur2_pink_2'; p05_load_sequ( p );
p.scan_name = 'ivw0033_Referenz_blau_5';p05_load_sequ( p );
p.scan_name = 'ivw0035_Struktur2_pink_3';p05_load_sequ( p );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Load force value: create figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.steps = [];
%p.out_path = '';'/gpfs/petra3/scratch/moosmanj/nextcloud';


p.raw_path = '/asap3/petra3/gpfs/p07/2021/data/11012199/raw';

p.steps = [];
p.out_path = '';
p.readhdf5 = 1;
p.adc2force = -19.92;

p.scan_name = 'ivw0015_referenzblau_1'; p05_load_force_values(p)
p.scan_name = 'ivw0015_referenzblau_2'; p05_load_force_values(p)
p.scan_name = 'ivw0015_referenzblau_2b'; p05_load_force_values(p)
p.scan_name = 'ivw0015_referenzblau_2c'; p05_load_force_values(p)
p.scan_name = 'ivw0016_referenzblau_3'; p05_load_force_values(p)
p.scan_name = 'ivw0017_Struktur1_gruen_1'; p05_load_force_values(p)
p.scan_name = 'ivw0018_Struktur1_gruen_1b'; p05_load_force_values(p)
p.scan_name = 'ivw0020_Struktur1_gruen_2'; p05_load_force_values(p)
p.scan_name = 'ivw0021_Struktur1_gruen_2b'; p05_load_force_values(p)
p.scan_name = 'ivw0022_Struktur1_gruen_2c'; p05_load_force_values(p)
p.scan_name = 'ivw0023_Struktur1_gruen_3'; p05_load_force_values(p)
p.scan_name = 'ivw0025_Struktur1_gruen_3'; p05_load_force_values(p)
p.scan_name = 'ivw0026_Struktur2_pink_1'; p05_load_force_values(p)
p.scan_name = 'ivw0027_Struktur2_pink_1b'; p05_load_force_values(p)
p.scan_name = 'ivw0027_Struktur2_pink_2'; p05_load_force_values(p)
p.scan_name = 'ivw0028_Struktur2_pink_1c'; p05_load_force_values(p)
p.scan_name = 'ivw0028_Struktur2_pink_1c'; p05_load_force_values(p)
p.scan_name = 'ivw0028_Struktur2_pink_2c'; p05_load_force_values(p)
p.scan_name = 'ivw0029_Struktur2_pink_2c_aperture05'; p05_load_force_values(p)
p.scan_name = 'ivw0030_Struktur2_pink_2c_aperture07'; p05_load_force_values(p)
p.scan_name = 'ivw0031_Struktur2_pink_2c_8001proj'; p05_load_force_values(p)
p.scan_name = 'ivw0032_Referenz_blau_4'; p05_load_force_values(p)
p.scan_name = 'ivw0033_Referenz_blau_5'; p05_load_force_values(p)
p.scan_name = 'ivw0035_Struktur2_pink_3'; p05_load_force_values(p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%