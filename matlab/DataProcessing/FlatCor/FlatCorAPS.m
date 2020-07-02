function FlatCorAPS(ParentPath,DataSet,TomoSetsToProcess,Cropping,HotPixThres_DarkFlatData,NumOfDarks,NumOfFlats,DarksBeforeFlats)
% Data preprocessing of the data from the 'Life-Cell Imaging' experiment
% taken at 2-BM@APS. Processing includes hot-pixel filtering darks, flats,
% and projections, dark-field correction, removal of stripe-like large-scale
% modulations in flats fields and projections which are due to the
% monochromator and are time-dependent, flat field correction of the such
% normatlized projections by the normalized flat fields

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters.
if nargin < 1
    % Assumed folder structure: 'ParentPath/data/DataSet'. Images of the
    % tomograms are found in 'DataSet'.
    ParentPath  = '/mnt/tomoraid3/tomo/APS_2BM_LifeCellImaging_GUP28266';
end
if nargin < 2
    % Name of input folder where images of the tomograms are found
    DataSet = 'wildtype_05min_deadtime_05tomo_stage16p0_upwards_620mm_050ms_30p0keV';
    %'/mnt/tomoraid3/tomo/APS_2BM_LifeCellImaging_GUP28266/data/wild_type_tomo_stage11p0_620mm_010ms_30p0keV';
end
if nargin < 3
    % Tomographic sets which should be processed. Scalar (0=default), or
    % 1x?-vector of numbers of the sets to process. Default processes all
    % sets found in 'DataSet'. 
    TomoSetsToProcess = 0;
end
if nargin < 4
    % scalar: 0=default cropping to [1008 1024], -1=check for croprange
    % file, 2x2-vector= [[HorFirst HorLast]; [VerFirst VerLast]]
    Cropping = -1;
end
if nargin < 5
    HotPixThres_DarkFlatData = [0.01 0.01 0.01];
end
if nargin < 6
    NumOfDarks = 60;
end
if nargin < 7
    NumOfFlats = 60;
end
if nargin < 8 
    DarksBeforeFlats = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Body.
tic;
DataPrefix = 'proj_';
DataFormat = 'tif';
NumOfProjs = 1200;
ProjOffset = 50;
NumOfImPerScan = 2*ProjOffset+NumOfProjs+NumOfDarks+NumOfFlats;
% Hot-pixel filter thresholds.
HotPixThresDark = HotPixThres_DarkFlatData(1);
HotPixThresFlat = HotPixThres_DarkFlatData(2);
HotPixThresData = HotPixThres_DarkFlatData(3);
%% input path
if ParentPath(end) ~= '/'
    ParentPath = [ParentPath '/'];
end
if DataSet(end) ~= '/'
    DataSet = [DataSet '/'];
end
InputPath     = [ParentPath 'data/' DataSet];
IntParentPath       = [ParentPath 'int/' DataSet];
OutputPathPre = [IntParentPath 'preprocessing/'];
%% Read folder and data structure.
DataStruct = dir([InputPath DataPrefix '*.' DataFormat]);
NumOfIm    = numel(DataStruct);
NumTomosFound = NumOfIm/(NumOfImPerScan);
% Print info.
fprintf('\nDATA SET TO PROCESS: %s\n',DataSet)
fprintf('INPUT PATH: %s\n',InputPath)
fprintf('TOTAL NUMBER OF IMAGES FOUND: %u\n',NumOfIm)
fprintf('CORRESPONDING NUMBER OF TOMOGRAMS: %u\n',NumTomosFound)
fprintf('TOMOGRAMS TO PROCESS:')
if TomoSetsToProcess == 0
    fprintf(' ALL\n')
else
    fprintf(' %s\n',mat2str(TomoSetsToProcess))
end
fprintf('SAVING ADDITIONAL IMAGES IN: %s',OutputPathPre)
if mod(NumTomosFound,1)~=0
    fprintf('\nIMAGES MISSING? NUMBER OF IMAGES FOUND IS NOT A MULTIPLE OF THE NUMBER OF IMAGES PER SCAN!\n')
end
% Read first image to get dimensions.
DataFileName = [InputPath DataStruct(1).name];
im           = double(imread(DataFileName,DataFormat));
% Data dimensions and optional cropping.
[dim1 dim2] = size(im);
dim1cen     = floor(dim1/2);
dim2cen     = floor(dim2/2);
%% Cropping.
if Cropping < 0
    Cropping = ReadCroprangeFile(IntParentPath);
end
if Cropping(1) == 0
    % Horizontal cropping.
    newYdim   = 1024;
    newYrange = dim2cen-newYdim/2+(1:newYdim);
    % Vertical cropping.
    newXdim   = 1008;
    newXrange = dim1cen-newXdim/2+(1:newXdim);
else
    % Horizontal cropping.
    newYrange = Cropping(1,1):Cropping(1,2);
    newYdim   = numel(newYrange);
    % Vertical cropping.
    newXrange = Cropping(2,1):Cropping(2,2);
    newXdim   = numel(newXrange);
end
if (newYdim>1) || (newXdim>1)
    fprintf('\nCROPPING IMAGES from [ver hor] = [%u %u] to [%u %u]!\n',dim1,dim2,newXdim,newYdim)
end
%% Loop over tomograms
if TomoSetsToProcess == 0
    TomoSetsToProcess = 1:NumTomosFound;
end
ttotal = toc;%disp(ttotal)
imcounter = 0;
for tomoNum = TomoSetsToProcess
    % Check output paths.
    OutputPath    = [ParentPath 'int/' DataSet sprintf('tomo%02u/',tomoNum)];
    fprintf('PROCESSING TOMO SET NUMBER: %u\n',tomoNum)
    fprintf('Saving processed data in: %s',OutputPath)
    % Creat output folder if not existing.
    if ~exist(OutputPath,'dir')
        mkdir(OutputPath);
    end
    if ~exist(OutputPathPre,'dir')
        mkdir(OutputPathPre);
    end
    %% Dark fields: reading, filtering, median.
    fprintf('\nReading dark fields:\n')
    tic;
    for nn = NumOfDarks:-1:1
        imIndex = (tomoNum-1)*NumOfImPerScan+NumOfProjs+2*ProjOffset+(1-DarksBeforeFlats)*NumOfFlats+nn;
        fprintf('%6u',imIndex)
        if mod(1+NumOfDarks-nn,30) == 0
            fprintf('\n')
        end
        DataFileName = [InputPath DataStruct(imIndex).name];
        im = double(imread(DataFileName,DataFormat));
        if (newYdim>1) || (newXdim>1)
            im = im(newXrange,newYrange);
        end
        im = FilterHotPixel(im,HotPixThresDark,0);
        imStack(:,:,nn) = im;
    end
    t = toc;ttotal = ttotal+t;
    [dimx dimy dimz] = size(imStack);
    if mod(1+NumOfDarks-nn,30) ~= 0
        fprintf('\n')
    end
    fprintf('Read %2u dark fields of dimension [%u %u] in %gs. Time to compute median: ',dimz,dimx,dimy,t)
    % Compute median of dark fields.
    tic;
    darkMed = median(imStack,3);
    t = toc; ttotal = ttotal+t;
    fprintf('%.2gs\n',t)
    edfwrite(sprintf('%sdark_fromMedianOfDarks_tomo%02u.edf',OutputPathPre,tomoNum),darkMed','float32');
    % Matrix for horizontal line integration (summation along y).
    tic;
    m = ones(dimy,1);
    %% Flat fields: reading, filtering, projecting, median.
    fprintf('Reading flat fields:\n')
    for nn = NumOfFlats:-1:1
        imIndex = (tomoNum-1)*NumOfImPerScan+NumOfProjs+2*ProjOffset+DarksBeforeFlats*NumOfDarks+nn;
        fprintf('%6u',imIndex)
        if mod(1+NumOfFlats-nn,30) == 0
            fprintf('\n')
        end
        DataFileName = [InputPath DataStruct(imIndex).name];
        im = double(imread(DataFileName,DataFormat));
        if (newYdim>1) || (newXdim>1)
            im = im(newXrange,newYrange);
        end
        im = FilterHotPixel(im-darkMed,HotPixThresFlat,0);
        %flat(:,:,nn) = im;
        flatLine(:,nn) = im*m/dimy;
        imStack(:,:,nn) = im./repmat(flatLine(:,nn),[1 dimy]);
    end
    t = toc;
    [dimx dimy dimz] = size(imStack);
    if mod(1+NumOfFlats-nn,30) ~= 0
        fprintf('\n')
    end
    fprintf('Read %2u flat fields of dimension [%u %u] in %gs. Time to compute median: ',dimz,dimx,dimy,t)
    % Compute median of dark fields.
    tic;
    flatMed = median(imStack,3);
    clear imStack
    t = toc; ttotal = ttotal+t;
    fprintf('%.2gs\n',t)
    edfwrite(sprintf('%sflat_fromMedianOfBgCorrectedFlats_tomo%02u.edf',OutputPathPre,tomoNum),flatMed','float32');
    edfwrite(sprintf('%sLineProjectedFlats_tomo%02u.edf',OutputPathPre,tomoNum),flatLine','float32');
    %% Read projections: reading, filtering, projecting, correlating.
    fprintf('Reading projections:\n')
    tic;
    for nn = NumOfProjs:-1:1
        imcounter = imcounter + 1;
        imIndex = (tomoNum-1)*NumOfImPerScan+ProjOffset+nn;
        fprintf('%6u',imIndex)
        if mod(1+NumOfProjs-nn,30) == 0
            fprintf('\n')
        end
        DataFileName = [InputPath DataStruct(imIndex).name];
        im = double(imread(DataFileName,DataFormat));
        if (newYdim>1) || (newXdim>1)
            im = im(newXrange,newYrange);
        end
        im = im - darkMed;
        %im = FilterHotPixel(im-darkMed,HotPixThresData,0);
        %data(:,:,nn) = im;
        dataLine(:,nn) = im*m/dimy;
        im = FilterHotPixel(im./repmat(dataLine(:,nn),[1 dimy])./flatMed,HotPixThresData,0);
        edfwrite(sprintf('%sint_%04u.edf',OutputPath,nn),im','float32');
    end
    t = toc;ttotal = ttotal+t;tic
    [dimx dimy] = size(im);
    if mod(1+NumOfProjs-nn,30) ~= 0
        fprintf('\n')
    end
    fprintf('Read %2u projections of dimension [%u %u] in %gs.\n',size(dataLine,2),dimx,dimy,t)
    edfwrite(sprintf('%sLineProjectedData_tomo%02u.edf',OutputPathPre,tomoNum),dataLine','float32');
    % Background shift of flat fields.
    for nn = NumOfFlats:-1:1
        vertShift(nn) = LineCorrelation(flatLine(:,1),flatLine(:,nn),0);
    end
    edfwrite(sprintf('%sVertShiftOfBgOfFlatsWrtFirstFlat_tomo%02u.edf',OutputPathPre,tomoNum),vertShift','float32');
    % Background shift of projections.
    for nn = NumOfProjs:-1:1
        vertShift(nn) = LineCorrelation(dataLine(:,1),dataLine(:,nn),0);
    end
    edfwrite(sprintf('%sVertShiftOfBgOfProjsWrtFirstProj_tomo%02u.edf',OutputPathPre,tomoNum),vertShift','float32');
    ttotal = ttotal+toc;
end
fprintf('\nA total of %u images processed in %g min, %g h, or %g s/images.\n',imcounter,ttotal/60,ttotal/3600,ttotal/imcounter)
