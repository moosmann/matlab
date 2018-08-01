% HNEE Experiment on wood with S. Lautner

%% Tomo reco loop
%/asap3/petra3/gpfs/common/p05/jm/matlab/experiments/p05/data/moosmanj/p05_reco_loop_hnee_tension_wood_000.m
edit p05_reco_loop_hnee_tension_wood_000


%% Load sequence processing
% volume slice mismatch
%[vol, vol_reg] = p05_load_sequ('/asap3/petra3/gpfs/p05/2018/data/11004450/processed','hnee23_pappel_tensionWood','reco_phase/tie_regPar1p50/float_rawBin2','x',3,0.0001, 0);

proc_path = '/asap3/petra3/gpfs/p05/2018/data/11004450/processed';
reco_subfolder = 'reco_phase/tie_regPar1p50/float_rawBin2';

scan = 'hnee18_pappel_tensionWood';
'hnee19_pappel_oppositeWood';
'hnee20_pappel_tensionWood';
'hnee21_pappel_oppositeWood';
'hnee23_pappel_tensionWood';
regdir = 'x';
steps = 2;[];
register = 1;
outlier_thresh = 0.0001;
auto_roi_center = 1;
[vol, vol_reg] = p05_load_sequ(proc_path, scan, reco_subfolder, regdir, steps, outlier_thresh, register, auto_roi_center);
