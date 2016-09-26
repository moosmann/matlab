function PipelineAPS(ParentPath)
%Start reconstruction pipeline of life-cell imaging data taken at 2-BM@APS.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Default arguments.
if nargin < 1
    % Assumed folder structure: 'ParentPath/data/DataSet'.
    ParentPath    = '/mnt/tomoraid-LSDF/tomo/APS_2BM_LifeCellImaging_GUP28266/lifecellimaging/';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prompt for input: data sets, tomo runs, preprocessing, phase retrieval, par-file creation, pyhst loop
DataSetStruct       = ChooseDatSet([ParentPath 'data/']);
NumDataSets         = numel(DataSetStruct);
TomoSetsToProcess   = input(sprintf('NUMBERS OF TOMOGRAMS TO RECONSTRUCT, ROW VECTOR, 0=ALL,DEFAULT: '));
if isempty(TomoSetsToProcess)
    TomoSetsToProcess = 0;
end
DoDataPreprocessing = input(sprintf('DATA PREPROCESSING, 1=YES, 0=NO,DEFAULT: '));
if isempty(DoDataPreprocessing)
    DoDataPreprocessing = 0;
end
DoPhaseRetrieval    = input(sprintf('PHASE RETRIEVAL, 1=YES, 0=NO,DEFAULT: '));
if isempty(DoPhaseRetrieval)
    DoPhaseRetrieval = 0;
end
DoCreatParFiles     = input(sprintf('CREATE PAR FILES, 1=YES, 0=NO,DEFAULT: '));
if isempty(DoCreatParFiles)
    DoCreatParFiles = 0;
end
if DoCreatParFiles == 1
    NumTomosForRotAxis = input(sprintf('TOMOGRAMS TO USE FOR ROTATION AXIS, 0=ALL,DEFAULT: '));
    if isempty(NumTomosForRotAxis)
        NumTomosForRotAxis = 0;
    end
    if NumTomosForRotAxis == 0
        NumTomosForRotAxis = TomoSetsToProcess;
    end
end
if isempty(DoCreatParFiles)
    DoCreatParFiles = 0;
end
DoPyhstLoop         = input(sprintf('LOOP PyHST OVER PAR FILES, 1=YES, 0=NO,DEFAULT: '));
if DoPyhstLoop == 1
    ParFilePrefix = input(sprintf('PART OF PAR FILE NAMES, STRING, ''''=ALL,''TIE''=DEFAULT: '));
    if isempty(ParFilePrefix)
        ParFilePrefix = 'TIE';
    end
else
    ParFilePrefix = 'TIE';
end
if isempty(DoPyhstLoop)
    DoPyhstLoop = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop over Data Sets.
fprintf('\nSTART RECONSTRUCTION PIPELINE')
for nn = 1:NumDataSets
    DataSet = DataSetStruct(nn).name;
    if DoDataPreprocessing == 1
        %% Data preprocessing: flat- and dark field correction, hot-pixel filtering.
        % Crop image before processing: 0: [1008 1024], or 2x2-vector
        Cropping = -2;
        filtmeth = 'LineSectionMed';
        medFiltLenHor = 63;
        HotPixThres_DarkFlatData = [0.001 0.001 0.001];
        FlatCorLCI(ParentPath,DataSet,TomoSetsToProcess,Cropping,filtmeth,medFiltLenHor,HotPixThres_DarkFlatData)
        %FlatCorAPS(ParentPath,DataSet,TomoSetsToProcess,DefaultCropping,HotPixThres_DarkFlatData)
    end
    %TomoSetsToProcess = 2:10;
    if DoPhaseRetrieval == 1
        %% Phase retrieval.
        alphaCTF_alphaTIE = 2.5;
        evalTIElo    = 1;
        evalTIEpnlo  = 0;
        evalCTF      = 0;
        BinaryFilterThreshold = 0;
        EnergyDistancePixelsize = [30 0.620 2.2e-6];
        InputFolderPrefix = 'tomo';
        PreProcessingFolder = 'int_filtLineSectionMedWidthH063V001_noCropping_filtSino';
        LCIFilterWidth_CircleCross = [60 00];
        %RecoLoopForAPS(alphaCTF_alphaTIE,evalTIElo,evalTIEpnlo,evalCTF,BinaryFilterThreshold,TomoSetsToProcess,EnergyDistancePixelsize,ParentPath,DataSet)
        RecoLoopLCI(alphaCTF_alphaTIE,evalTIElo,evalTIEpnlo,evalCTF,BinaryFilterThreshold,TomoSetsToProcess,EnergyDistancePixelsize,ParentPath,DataSet,InputFolderPrefix,PreProcessingFolder,LCIFilterWidth_CircleCross)
    end
    %% Make par files for PyHST.
    %TomoSetsToProcess = 1:10;
    if DoCreatParFiles == 1
        StartEndVoxels = -1;
        NumOfFirstAndLastProjection = [1 1200];
        EffectivePixelSize = 2.2;
        AngleBetweenProjections = 180/1200;
        InputFolderPrefix = 'tomo';
        PreProcessingFolder = 'int_filtLineSectionMedWidthH063V001_noCropping_filtSino';
        MakeParFileLCI(StartEndVoxels,NumOfFirstAndLastProjection,EffectivePixelSize,AngleBetweenProjections,ParentPath,DataSet,NumTomosForRotAxis,InputFolderPrefix,PreProcessingFolder)
    end
    %% Start tomographic reconstruciton using PyHST.
    if DoPyhstLoop == 1
        VolPath = [ParentPath 'vol/' DataSet];
        Pyhst(VolPath,ParFilePrefix,0)
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('RECONSTRUCTION PIPELINE FINISHED.\n')

