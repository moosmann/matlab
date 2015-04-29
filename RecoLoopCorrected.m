clear all
%% Data set: wildtype_30keV_10min_deadtime_10tomo_stage11p0_upwards_620mm_050ms
%% Set parameters
alphaCTF_alphaTIE = 3;
evalTIElo    = 1;
evalTIEpnlo  = 0;
evalCTF      = 0;
BinaryFilterThreshold = 0;
EnergyPixelsizeDistance = [30 2.2e-6 0.620];
InputPath = '/mnt/tomoraid-LSDF/tomo/APS_2BM_LifeCellImaging_GUP28266/lifecellimaging/preprocess/corrected';
FileNamePrefix = 'proj';
OutputPath = '/mnt/tomoraid-LSDF/tomo/APS_2BM_LifeCellImaging_GUP28266/lifecellimaging/preprocess/corrected_phase';
doSino = 0;
%% Start phase retrieval loop
RecoLoopThis(alphaCTF_alphaTIE,evalTIElo,evalTIEpnlo,evalCTF,BinaryFilterThreshold,EnergyPixelsizeDistance,InputPath,FileNamePrefix,OutputPath,doSino);