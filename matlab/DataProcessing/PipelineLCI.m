function PipelineLCI()
%Start reconstruction pipeline of life-cell imaging data taken at 2-BM@APS.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ParentPath    = '/mnt/tomoraid-LSDF/tomo/APS_2BM_LifeCellImaging_GUP28266/lifecellimaging/';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prompt for input: data sets, tomo runs, preprocessing, phase retrieval, par-file creation, pyhst loop
DoDataPreprocessing = input(sprintf('DATA PREPROCESSING, 0=NO,DEFAULT, 1=GAUSSIAN, 2=BUTTERWORTH, 3=LINE: '));
if isempty(DoDataPreprocessing)
    DoDataPreprocessing = 0;
end
if DoDataPreprocessing > 0
    if DoDataPreprocessing == 1
        GaussianWidth = input(sprintf('  Gaussian width for Fourier space filter: '));
    elseif DoDataPreprocessing == 2
        ButterworthParameter = input(sprintf('  Butterworth filter parameters: [RelativeCutoff FilterOrder] = '));
    elseif DoDataPreprocessing == 3
        filtmeth = input(sprintf('  Filter method: 1=LineSectionMed, 2=LineCompleteMean 3=LineCompleteMed: '));
        switch filtmeth
            case 1
                filtmeth = 'LineSectionMed';
            case 2
                filtmeth = 'LineCompleteMean';
            case 3
                filtmeth = 'LineCompleteMed';
        end
        if strcmp(filtmeth,'LineSectionMed')
            medFiltLenHor = input(sprintf('  Horizontal size of median filter: '));
        end
    end
    doDarkFieldCorrection = input(sprintf('  Apply dark field correction, 0=no, 1=yes: '));
end
DataSetStruct       = ChooseDatSet([ParentPath '/int']);
NumDataSets         = numel(DataSetStruct);
TomoSetsToProcess   = input(sprintf('NUMBERS OF TOMOGRAMS TO RECONSTRUCT, ROW VECTOR, 0=ALL,DEFAULT: '));
if isempty(TomoSetsToProcess)
    TomoSetsToProcess = 0;
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
    TomosForRotAxis = input(sprintf('TOMOGRAMS TO USE FOR ROTATION AXIS, 0=ALL,DEFAULT: '));
    if isempty(TomosForRotAxis)
        TomosForRotAxis = 0;
    end
    if TomosForRotAxis == 0
        TomosForRotAxis = TomoSetsToProcess;
    end
end
if isempty(DoCreatParFiles)
    DoCreatParFiles = 0;
end
DoPyhstLoop         = input(sprintf('LOOP PyHST OVER PAR FILES, 1=YES, 0=NO,DEFAULT: '));
if DoPyhstLoop == 1
    ParFilePrefix = input(sprintf('  Part of par file names, string, ''''=ALL,''TIE''=DEFAULT: '));
    if isempty(ParFilePrefix)
        ParFilePrefix = 'TIE';
    end
else
    ParFilePrefix = 'TIE';
end
if isempty(DoPyhstLoop)
    DoPyhstLoop = 0;
end
fprintf('\nSTART RECONSTRUCTION PIPELINE')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nn = 1:NumDataSets
    DataSet = DataSetStruct(nn).name;
    PreProcessingFolder = '';
    if strcmp(DataSet,'wildtype_30keV_10min_deadtime_20tomo_stage14p0_upwards_620mm_025ms') || strcmp(DataSet,'wildtype_30keV_10min_deadtime_25tomo_stage11p0_upwards_620mm_015ms')
        PreProcessingFolder = 'int_filtLineSectionMedWidthH063V001_noCropping';        
        PreProcessingFolder = [PreProcessingFolder '_filtSino'];
    end
    disp(PreProcessingFolder);
    if DoDataPreprocessing > 0
        %% Data preprocessing: flat- and dark field correction, hot-pixel filtering.
        % Crop image before processing: 0: [1008 1024], or 2x2-vector
        Cropping = -2;
        HotPixThres_DarkFlatData = [1 1 1];%[0.04 0.04 0.04];%[0.2 0.2 0.2];%
        if DoDataPreprocessing == 1
            FlatCorLCI2(ParentPath,DataSet,TomoSetsToProcess,Cropping,GaussianWidth,HotPixThres_DarkFlatData,doDarkFieldCorrection)
        elseif DoDataPreprocessing == 2
            FlatCorLCI3(ParentPath,DataSet,TomoSetsToProcess,Cropping,ButterworthParameter,HotPixThres_DarkFlatData,doDarkFieldCorrection)
        elseif DoDataPreprocessing == 3
            FlatCorLCI(ParentPath,DataSet,TomoSetsToProcess,Cropping,filtmeth,medFiltLenHor,HotPixThres_DarkFlatData,doDarkFieldCorrection)
        end
    end
    % Phase retrieval
    if DoPhaseRetrieval == 1
        doRingFilter = 0;
        %RecoLoopLCI(alphaCTF_alphaTIE,evalTIElo,evalTIEpnlo,evalCTF,BinaryFilterThreshold,TomoSetsToProcess,EnergyDistancePixelsize,ParentPath,DataSet,InputFolderPrefix,PreProcessingFolder);
        RecoLoopLCIstack(TomoSetsToProcess,doRingFilter,'ParentPath',ParentPath,'DataSet',DataSet,'PreProcessingFolder',PreProcessingFolder);
    end
    %% Make par files for PyHST
    if DoCreatParFiles == 1
        StartEndVoxels = -1;
        NumOfFirstAndLastProjection = [1 1200];
        EffectivePixelSize = 2.2;
        AngleBetweenProjections = 180/1200;
        MakeParFileLCI(StartEndVoxels,NumOfFirstAndLastProjection,EffectivePixelSize,AngleBetweenProjections,ParentPath,DataSet,TomosForRotAxis,PreProcessingFolder)
    end
    %% Start tomographic reconstruction using PyHST.
    if DoPyhstLoop == 1
        VolPath = [ParentPath '/vol/' PreProcessingFolder '/' DataSet];
        Pyhst(VolPath,ParFilePrefix,0)
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('RECONSTRUCTION PIPELINE FINISHED.\n')