% HNEE Experiment on wood with S. Lautner
% 2018 April 20, 11004450

%% Tomo reco loop
%/asap3/petra3/gpfs/common/p05/jm/matlab/experiments/p05/data/moosmanj/p05_reco_loop_hnee_tension_wood_000.m
edit p05_reco_loop_hnee_tension_wood_000
edit p05_reco_loop_hnee.m

%% Load sequence processing
p.proc_path = '/asap3/petra3/gpfs/p05/2018/data/11004450/processed';
p.reco_sub = 'reco_phase/tie_regPar1p50/float_rawBin2_recoBin2';
p.steps = [];
p.out_thresh = 0.0001;

% GERISSEN, 9 scans, 8 steps
p.scan_name = 'hnee18_pappel_tensionWood';
p.roi = [116 847 185 419 135 703];

% Kleber Probleme
p.scan_name = 'hnee20_pappel_tensionWood';
p.roi = [137 867 403 578 246 678];

% GERISSEN nach 5. scan, bei ca 35 N
p.scan_name = 'hnee21_pappel_oppositeWood';
p.roi = [79 927 393 557 255 680];

p.scan_name = 'hnee23_pappel_tensionWood';
p.roi = [144 820 613 735 252 680];

p.scan_name = 'hnee22_pappel_oppositeWood';
p.roi = [];

% GERISSEN, vor 8. scan
p.scan_name = 'hnee19_pappel_oppositeWood';
p.roi = [225 890 355 490 171 690];

% GERISSEN nach 5. scan, bei ca 35 N
p.scan_name = 'hnee21_pappel_oppositeWood';
p.roi = [79 927 393 557 255 680];

p.reco_sub = 'reco_phase/tie_regPar1p50/float_rawBin2';
p.steps = [5 7];
[vol, roi] = p05_load_sequ_2( p );

%% Load sequence force plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.raw_path = '/asap3/petra3/gpfs/p05/2018/data/11004450/raw';
p.steps = [];
p.out_path = '';
% 20 N
weight = 5.31;
voltage = 5.08 / 2;
gconst = 9.81;
p.adc2force = -weight * gconst / voltage;
p.readhdf5 = 1;

% GERISSEN, 9 scans, 8 steps
p.scan_name = 'hnee18_pappel_tensionWood';
%p05_load_force_values( p )

% GERISSEN, vor 8. scan %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Good one, broke during wait before last scan
p.scan_name = 'hnee19_pappel_oppositeWood';
%p05_load_force_values( p )

p.scan_name = 'hnee20_pappel_tensionWood';
p05_load_force_values( p )

% GERISSEN nach 5. scan, bei ca 35 N
p.scan_name = 'hnee21_pappel_oppositeWood';
%p05_load_force_values( p )

p.scan_name = 'hnee22_pappel_oppositeWood';
p05_load_force_values( p )

p.scan_name = 'hnee23_pappel_tensionWood';
p05_load_force_values( p )


