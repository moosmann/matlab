%% Directories, folder, parameters.
% Assumed folder structure: 'ParentPath/data/DataSet'.
ParentPath  = '/mnt/tomoraid-LSDF/tomo/APS_2BM_LifeCellImaging_GUP28266/lifecellimaging/';
%'/mnt/tomoraid3/tomo/APS_2BM_LifeCellImaging_GUP28266';
DataSet = 'wildtype_30keV_10min_deadtime_25tomo_stage11p0_upwards_620mm_015ms';
%'wildtype_30keV_10min_deadtime_25tomo_stage11p0_upwards_620mm_015ms';
%'wildtype_05min_deadtime_05tomo_stage16p0_upwards_620mm_050ms_30p0keV';
%'wild_type_tomo_stage11p0_620mm_010ms_30p0keV';
% 0 = all found.
TomoSetsToProcess = 23:25;
fprintf('\nSTART RECONSTRUCTION PIPELINE\n')
%% Data preprocessing: flat- and dark field correction, hot-pixel filtering.
% Crop image before processing: 0: [1008 1024], or 2x2-vector
DefaultCropping = -1;
HotPixThres_DarkFlatData = [0.01 0.01 0.01];
FlatCorAPS(ParentPath,DataSet,TomoSetsToProcess,DefaultCropping,HotPixThres_DarkFlatData)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TomoSetsToProcess = 1:22;
%% Phase retrieval.
alphaCTF_alphaTIE = 2.5;
evalTIElo    = 1;
evalTIEpnlo  = 0;
evalCTF      = 0;
BinaryFilterThreshold = 0;
EnergyDistancePixelsize = [30 0.620 2.2e-6];
RecoLoopForAPS(alphaCTF_alphaTIE,evalTIElo,evalTIEpnlo,evalCTF,BinaryFilterThreshold,TomoSetsToProcess,EnergyDistancePixelsize,ParentPath,DataSet)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make par files for PyHST.
%TomoSetsToProcess = 1:10;
StartEndVoxels = -1;
NumOfFirstAndLastProjection = [1 1200];
EffectivePixelSize = 2.2;
AngleBetweenProjections = 180/1200;
NumTomosForRotAxis = 21;
MakeParFileForAPS(StartEndVoxels,NumOfFirstAndLastProjection,EffectivePixelSize,AngleBetweenProjections,ParentPath,DataSet,NumTomosForRotAxis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start tomographic reconstruciton using PyHST.
VolPath = [ParentPath 'vol/' DataSet];
Pyhst(VolPath,'',0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%