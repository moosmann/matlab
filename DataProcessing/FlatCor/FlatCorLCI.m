function FlatCorLCI(ParentPath,DataSet,TomoSetsToProcess,Cropping,filtmeth,medFiltLenHor,HotPixThres_DarkFlatData,doDarkFieldCorrection,NumOfDarks,NumOfFlats,DarksBeforeFlats)
% Data preprocessing of the data from the 'Life-Cell Imaging' experiment
% taken at 2-BM@APS. Processing includes hot-pixel filtering darks, flats,
% and projections, dark-field correction, removal of stripe-like large-scale
% modulations in flats fields and projections which are due to the
% monochromator and are time-dependent, flat field correction of the such
% normatlized projections by the normalized flat fields
% Modfied to include line section filtering and for all data set on
% 29/08/2012

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters.
if nargin < 1
    % Assumed folder structure: 'ParentPath/data/DataSet'. Images of the
    % tomograms are found in 'DataSet'.
    ParentPath = '/mnt/tomoraid-LSDF/tomo/APS_2BM_LifeCellImaging_GUP28266/lifecellimaging/';
    %'/mnt/tomoraid-LSDF/tomo/APS_2BM_LifeCellImaging_GUP28266/lifecellimaging';
end
if nargin < 2
    % Name of input folder where images of the tomograms are found
    DataSet = 'wildtype_30keV_10min_deadtime_25tomo_stage11p0_upwards_620mm_015ms';
    %'/mnt/tomoraid3/tomo/APS_2BM_LifeCellImaging_GUP28266/data/wild_type_tomo_stage11p0_620mm_010ms_30p0keV';
end
if nargin < 3
    % Tomographic sets which should be processed. Scalar (0=default), or
    % 1x?-vector of numbers of the sets to process. Default processes all
    % sets found in 'DataSet'. 
    TomoSetsToProcess = 1;
end
if nargin < 4
    % scalar: 0=default cropping to [1008 1024], -1=check for croprange
    % file, 2x2-vector= [[HorFirst HorLast]; [VerFirst VerLast]]
    Cropping = -2;
end
if nargin < 5
    filtmeth = 'LineSectionMed';
end
if nargin < 6
    medFiltLenHor = 63;
end
if nargin < 7
    HotPixThres_DarkFlatData = [0.001 0.001 0.001];
end
if nargin < 8
    doDarkFieldCorrection = 1;
end
if nargin < 9
    NumOfDarks = 60;
end
if nargin < 10
    NumOfFlats = 60;
end
if nargin < 11
    DarksBeforeFlats = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if filtmeth == 0
    filtmeth = input(sprintf('  Type in as a string one of the following filter methods, LineSectionMed LineCompleteMean LineCompleteMed: '));
end
if strcmp(filtmeth,'LineSectionMed')
    if medFiltLenHor == 0
        medFiltLenHor = input(sprintf('  Type in size of median filter (in horizontal direction): '));
    end
end
medFiltLenVer = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Body.
tic;
DataPrefix = 'proj_';
DataFormat = 'tif';
NumOfProjs = 1200;
ProjOffset = 50;
NumOfImPerScan = 2*ProjOffset+NumOfProjs+NumOfDarks+NumOfFlats;
% Aterfact filtering by 2D median filtering.
x1 = 223;
y1 = 389;
x2 = 237;
y2 = 403;
medfiltwidth_x = ceil((x2-x1)/2+1);
medfiltwidth_y = ceil((y2-y1)/2+1);
% Hot-pixel filter thresholds.
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
if Cropping == 0
    cropmeth = '_defaultCropping';
elseif Cropping == -1
    cropmeth = '_cropped';
elseif Cropping == -2
    cropmeth = '_noCropping';
end
if doDarkFieldCorrection == 1
    dodark = '';
elseif doDarkFieldCorrection == 0
    dodark = '_noDarkFieldCorrection';
end
if isequal(filtmeth,'LineSectionMed')
    IntNamePostfix = ['filt' filtmeth sprintf('WidthH%03uV%03u',medFiltLenHor,medFiltLenVer) cropmeth dodark];
else
    IntNamePostfix = ['filt' filtmeth cropmeth dodark];
end
IntName = ['int_' IntNamePostfix];
OutputParentPath = '/mnt/tomoraid-LSDF/tomo/APS_2BM_LifeCellImaging_GUP28266/lifecellimaging/int/';
OutputPathPre = [OutputParentPath DataSet '/' IntName '/preprocessing/'];
%% Read folder and data structure.
tic;
DataStruct = dir([InputPath DataPrefix '*.' DataFormat]);
fprintf('Time to read data struct: %f s\n',toc);
NumOfIm    = numel(DataStruct);
NumTomosFound = NumOfIm/(NumOfImPerScan);
% Print info.
fprintf('DATA SET TO PROCESS: %s\n',DataSet)
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
im           = single(imread(DataFileName,DataFormat));
% Data dimensions and optional cropping.
[dim1,dim2] = size(im);
dim1cen     = floor(dim1/2);
dim2cen     = floor(dim2/2);
%% Cropping.
if Cropping == -1
    Cropping = ReadCroprangeFile(IntParentPath);
end
if Cropping(1) == 0
    % Horizontal cropping.
    newYdim   = 1024;
    newYrange = dim2cen-newYdim/2+(1:newYdim);
    % Vertical cropping.
    newXdim   = 1008;
    newXrange = dim1cen-newXdim/2+(1:newXdim);
elseif Cropping(1) == -2
    newXdim = 0;
    newYdim = 0;
else
    % Horizontal cropping.
    newYrange = Cropping(1,1):Cropping(1,2);
    newYdim   = numel(newYrange);
    % Vertical cropping.
    newXrange = Cropping(2,1):Cropping(2,2);
    newXdim   = numel(newXrange);
end
if (newYdim>1) || (newXdim>1)
    [~, dimy] = size(im(newXrange,newYrange));
    fprintf('\nCROPPING IMAGES from [ver hor] = [%u %u] to [%u %u]!\n',dim1,dim2,newXdim,newYdim)
else
    %dimx = dim1;
    dimy = dim2;
end
%% Loop over tomograms
if TomoSetsToProcess == 0
    TomoSetsToProcess = 1:NumTomosFound;
end
ttotal = toc;%disp(ttotal)
imcounter = 0;
for tomoNum = TomoSetsToProcess
    % Check output paths.
    OutputPath    = [OutputParentPath DataSet IntName '/' sprintf('tomo%02u/',tomoNum)];
    fprintf('\nPROCESSING TOMO SET NUMBER: %u\n',tomoNum)
    fprintf('Saving processed data in: %s',OutputPath)
    % Creat output folder if not existing.
    if ~exist(OutputPath,'dir')
        mkdir(OutputPath);
    end
    if ~exist(OutputPathPre,'dir')
        mkdir(OutputPathPre);
    end
    %% DARK FIELDS: reading, filtering, median.
    if doDarkFieldCorrection
        HotPixThresDark = HotPixThres_DarkFlatData(1);
        fprintf('\nReading dark fields:\n')
        tic;
        for nn = NumOfDarks:-1:1
            imIndex = (tomoNum-1)*NumOfImPerScan+NumOfProjs+2*ProjOffset+(1-DarksBeforeFlats)*NumOfFlats+nn;
            fprintf('%6u',imIndex)
            if mod(1+NumOfDarks-nn,30) == 0
                fprintf('\n')
            end
            DataFileName = [InputPath DataStruct(imIndex).name];
            im = imread(DataFileName,DataFormat);
            if (newYdim>1) || (newXdim>1)
                im = im(newXrange,newYrange);
            end
            im = FilterHotPixel(im,HotPixThresDark,0);
            imStack(:,:,nn) = im;
        end
        t = toc;ttotal = ttotal+t;
        [dimx, dimy, dimz] = size(imStack);
        if mod(1+NumOfDarks-nn,30) ~= 0
            fprintf('\n')
        end
        fprintf('Read %2u dark fields of dimension [%u %u] in %gs. Time to compute median: ',dimz,dimx,dimy,t)
        % Compute median of dark fields.
        tic;
        darkMed = single(median(imStack,3));
        clear imStack
        fprintf('%.2gs\n',t)
        edfwrite(sprintf('%sdark_fromMedianOfDarks_tomo%02u_%s.edf',OutputPathPre,tomoNum,IntNamePostfix),darkMed','float32');
        t = toc; ttotal = ttotal+t;
    end
    % Matrix for horizontal line integration (summation along y).
    tic;
    m = ones(dimy,1);
    %% FLAT FIELDs: reading, filtering, projecting, median.
    fprintf('Reading flat fields:\n')
    for nn = NumOfFlats:-1:1
        imIndex = (tomoNum-1)*NumOfImPerScan+NumOfProjs+2*ProjOffset+DarksBeforeFlats*NumOfDarks+nn;
        fprintf('%6u',imIndex)
        if mod(1+NumOfFlats-nn,30) == 0
            fprintf('\n')
        end
        DataFileName = [InputPath DataStruct(imIndex).name];
        im = single(imread(DataFileName,DataFormat));
        if (newYdim>1) || (newXdim>1)
            im = im(newXrange,newYrange);
        end
        im = FilterHotPixel(im,HotPixThresFlat,0);
        if exist('darkMed','var')
            im = im - darkMed;
        end
        switch filtmeth
            case 'LineCompleteMean'
                %flat(:,:,nn) = im;
                flatLine(:,nn)  = im*m/dimy;
                imStack(:,:,nn) = im./repmat(flatLine(:,nn),[1 dimy]);
            case 'LineCompleteMed'
                %flat(:,:,nn) = im;
                flatLine(:,nn)  = median(im,2);
                imStack(:,:,nn) = im./repmat(flatLine(:,nn),[1 dimy]);
            case 'LineSectionMed'
                imMedFilt = medfilt2(im,[medFiltLenVer medFiltLenHor],'symmetric');
                imStack(:,:,nn) = im./imMedFilt;
            case 'None'
                imStack(:,:,nn) = im;
        end
    end
    t = toc;
    [dimx, dimy, dimz] = size(imStack);
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
    edfwrite(sprintf('%sflat_fromMedianOfBgCorrectedFlats_tomo%02u_%s.edf',OutputPathPre,tomoNum,IntNamePostfix),flatMed','float32');
    if exist('flatLine','var')
        edfwrite(sprintf('%sLineProjectedFlats_tomo%02u_%s.edf',OutputPathPre,tomoNum,IntNamePostfix),flatLine','float32');
    end
    %% PROJECTIONs: reading, filtering, projecting, correlating.
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
        im = single(imread(DataFileName,DataFormat));
        if (newYdim>1) || (newXdim>1)
            im = im(newXrange,newYrange);
        end
        im = FilterHotPixel(im,HotPixThresData,0);
        if exist('darkMed','var')
            im = im - darkMed;
        end
        %im = FilterHotPixel(im-darkMed,HotPixThresData,0);
        switch filtmeth
            case 'LineCompleteMean'
                dataLine(:,nn) = im*m/dimy;
                im = im./repmat(dataLine(:,nn),[1 dimy])./flatMed;
            case 'LineCompleteMed'
                dataLine(:,nn) = median(im,2);
                im = im./repmat(dataLine(:,nn),[1 dimy])./flatMed;
            case 'LineSectionMed'
                im = im./medfilt2(im,[medFiltLenVer medFiltLenHor],'symmetric');
                im = im./flatMed;
            case 'None'
                im = im./flatMed;
        end
        % Filter artefact.
        im(x1:x2,y1:y2) = medfilt2(im(x1:x2,y1:y2),[medfiltwidth_x medfiltwidth_y],'symmetric');
        edfwrite(sprintf('%sint_%04u.edf',OutputPath,nn),im','float32');
    end
    t = toc;ttotal = ttotal+t;tic
    [dimx, dimy] = size(im);
    if mod(1+NumOfProjs-nn,30) ~= 0
        fprintf('\n')
    end
    if exist('dataLine','var')
        fprintf('Read %2u projections of dimension [%u %u] in %gs.\n',size(dataLine,2),dimx,dimy,t)
        edfwrite(sprintf('%sLineProjectedData_tomo%02u_%s.edf',OutputPathPre,tomoNum,IntNamePostfix),dataLine','float32');
        % Background shift of flat fields.
        for nn = NumOfFlats:-1:1
            vertShift(nn) = LineCorrelation(flatLine(:,1),flatLine(:,nn),0);
        end
        % Background shift of projections.
        for nn = NumOfProjs:-1:1
            vertShift(nn) = LineCorrelation(dataLine(:,1),dataLine(:,nn),0);
        end
        edfwrite(sprintf('%sVertShiftOfBgOfProjsWrtFirstProj_tomo%02u_%s.edf',OutputPathPre,tomoNum,IntNamePostfix),vertShift','float32');
    end
    ttotal = ttotal+toc;
end
fprintf('\nA total of %u images processed in %g min, %g h, or %g s/images.\n',imcounter,ttotal/60,ttotal/3600,ttotal/imcounter)
