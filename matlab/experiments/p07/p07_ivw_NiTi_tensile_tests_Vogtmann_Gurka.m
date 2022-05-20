% Kaierslautern Leibniz Insitute IVW
% Julia Vogtmann, Martin Gurka
% Nickel Titan wires
%
% First BT: energy 65 keV
% Second beamtime: reduce energy

edit p07_ivw_NiTi_tensile_tests_Vogtmann_Gurka.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SECOND BT 11012816 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

edit p07_ivw_20220304_11012816.m

%% Extract scan names to be stitched %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\n\nIDENTIFY SCANS TO BE STITCHED:' )
proc_path = sprintf('/asap3/petra3/gpfs/p07/2022/data/11012816/processed/');
d = dir( proc_path );
s = {};
sn = 1;
for n = 1:numel(d)
    name1 = d(n).name;
    if strcmp( name1(end), 'a')
        name2 = d(n+1).name;
        if strcmp( name2(end), 'b')
            name0 = name1(1:end-2);            
            s{sn} = name0;            
            fprintf( '\n %4u %s %s %s', sn, name1, name2, name0)
            sn = sn + 1;   
        end
    end
end

%% Stitch scans
fprintf( '\n\nSTITCH SCANS' )
for n = 26%1:numel(s) 
    name = s{n};
    scan_path = sprintf( '%s%s', proc_path, name);
    scan_subfolder = 'reco';
    reco_subfolder = 'float_rawBin3';
    crop = 1;
    save_stitched_volume = 1;    
    stitched_volume_path = '';
    
    full_path_a = sprintf( '%s_a/%s/%s', scan_path, scan_subfolder, reco_subfolder );
    full_path_b = sprintf( '%s_b/%s/%s', scan_path, scan_subfolder, reco_subfolder );
    
    fprintf( '\n\n %4u %s %u %u', n, scan_path, exist( full_path_a, 'dir' ), exist( full_path_b, 'dir' ) )
    
    %[s, vol] = stitch_volumes( scan_path, scan_subfolder, reco_subfolder, crop, save_stitched_volume, stitched_volume_path );
    stitch_volumes( scan_path, scan_subfolder, reco_subfolder, crop, save_stitched_volume, stitched_volume_path );
end

fprintf( '\n' )

%% Load force value: create figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Extract scan names  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\n\nIDENTIFY SCANS and EXTRACT LOAD VALUES:' )
proc_path = sprintf('/asap3/petra3/gpfs/p07/2022/data/11012816/raw/*_000_setforce');
d = dir( proc_path );
p.steps = [];
%p.out_path = '';
p.raw_path = '/asap3/petra3/gpfs/p07/2022/data/11012816/raw';
p.out_path = '';
p.readhdf5 = 1;
p.adc2force = [];

for n=1:numel(d)
    p.scan_name = d(n).name(1:end-13);
    fprintf('\n %2u: %s', n, p.scan_name )
    load_force_values(p)
end
fprintf( '\nLOOP FINISHED' )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIRST BT 11012199 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

edit p07_ivw_20210701_11012199.m

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
load_sequ( p );

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
vol = load_sequ( p );
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
load_sequ( p );

%p.scan_name = 'ivw0017_Struktur1_gruen_1'; load_sequ( p );
p.scan_name = 'ivw0020_Struktur1_gruen_2'; load_sequ( p );
p.scan_name = 'ivw0021_Struktur1_gruen_2b';load_sequ( p );
p.scan_name = 'ivw0025_Struktur1_gruen_3';load_sequ( p );
p.scan_name = 'ivw0027_Struktur2_pink_1b'; load_sequ( p );
p.scan_name = 'ivw0027_Struktur2_pink_2'; load_sequ( p );
p.scan_name = 'ivw0033_Referenz_blau_5';load_sequ( p );
p.scan_name = 'ivw0035_Struktur2_pink_3';load_sequ( p );
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

p.scan_name = 'ivw0015_referenzblau_1'; load_force_values(p)
p.scan_name = 'ivw0015_referenzblau_2'; load_force_values(p)
p.scan_name = 'ivw0015_referenzblau_2b'; load_force_values(p)
p.scan_name = 'ivw0015_referenzblau_2c'; load_force_values(p)
p.scan_name = 'ivw0016_referenzblau_3'; load_force_values(p)
p.scan_name = 'ivw0017_Struktur1_gruen_1'; load_force_values(p)
p.scan_name = 'ivw0018_Struktur1_gruen_1b'; load_force_values(p)
p.scan_name = 'ivw0020_Struktur1_gruen_2'; load_force_values(p)
p.scan_name = 'ivw0021_Struktur1_gruen_2b'; load_force_values(p)
p.scan_name = 'ivw0022_Struktur1_gruen_2c'; load_force_values(p)
p.scan_name = 'ivw0023_Struktur1_gruen_3'; load_force_values(p)
p.scan_name = 'ivw0025_Struktur1_gruen_3'; load_force_values(p)
p.scan_name = 'ivw0026_Struktur2_pink_1'; load_force_values(p)
p.scan_name = 'ivw0027_Struktur2_pink_1b'; load_force_values(p)
p.scan_name = 'ivw0027_Struktur2_pink_2'; load_force_values(p)
p.scan_name = 'ivw0028_Struktur2_pink_1c'; load_force_values(p)
p.scan_name = 'ivw0028_Struktur2_pink_1c'; load_force_values(p)
p.scan_name = 'ivw0028_Struktur2_pink_2c'; load_force_values(p)
p.scan_name = 'ivw0029_Struktur2_pink_2c_aperture05'; load_force_values(p)
p.scan_name = 'ivw0030_Struktur2_pink_2c_aperture07'; load_force_values(p)
p.scan_name = 'ivw0031_Struktur2_pink_2c_8001proj'; load_force_values(p)
p.scan_name = 'ivw0032_Referenz_blau_4'; load_force_values(p)
p.scan_name = 'ivw0033_Referenz_blau_5'; load_force_values(p)
p.scan_name = 'ivw0035_Struktur2_pink_3'; load_force_values(p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%