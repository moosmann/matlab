function DataProc32IDstack201307tomo(varargin)
% data processing for the experiment In vivo imaging, GUP at beamline
% 32ID-C at APS, ANL, Chicago, Illinois, USA, 2013-07-28 to 2013-08-02.
% Allows for tomographic reconstruction. Hower, due to Matlab's slow
% interpolation function, the use of PyHST is much faster. For data
% processing including gathering of hot pixel statistics, see
% DataProc32IDstack201307tomoWithStatistics. FBP option deleted, see
% DataProc32IDstack2013tomoWithFBP.
%
%Written by Julian Moosmann, 2013-09-19

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
%% PARSE INPUT ARGUMENT 'VARARGIN'
ParseVarargin(varargin)
%% DEFAULT ARGUMENTS AND PARAMETERS
if ~exist('TomoToProcess','var')
    % Index of tomogram as in filename. It's not the
    % absolute number of tomograms!
    TomoToProcess = 1;
end
if ~exist('ParentPath','var')
    ParentPath = '/mnt/tomoraid-LSDF/tomo/APS_32ID-C_InVivoImaging_GUP34879_2013-07-30/';
end
if ~exist('DataSet','var')
    %DataSet = 'Lamprey_atlas/Aug02_01-01_Lamprey_stage21p0_22p0keV_0700mm_50ms_1800proj';
    %DataSet = 'Lamprey_inVivo/August01_07-10_Lamprey_stage10somites__30p0keV_0700mm_15ms_500proj_scantime20s_deadtime08min';
    %DataSet = 'Zebrafish_inVivo/July31_20-00_Zebrafish_stage07p0hpf_30p0keV_0700mm_15ms_500proj_scantime20s_deadtime08min';
    %DataSet = 'Zebrafish_atlas/August02_02-30_Zebrafish_stage20p0hpf_22p0keV_0700mm_50ms_1800proj';
    %DataSet = 'Xenopus_inVivo/Jul28_14-50_urea_stage17p0_30p0keV_0400mm_20ms_0500proj_scantime20s_deadtime8min';
    %DataSet = 'Xenopus_inVivo/Jul28_20-25_urea_stage21p0_30p0keV_0400mm_10ms_0500proj_scantime20s_deadtime8min';
    %DataSet = 'Xenopus_inVivo/Jul29_05-15_urea_stage25p0_30p0keV_0400mm_10ms_0500proj_scantime20s_deadtime8min';
    %DataSet = 'Xenopus_inVivo/Jul29_15-10_urea_stage27p0_30p0keV_0700mm_15ms_0500proj_scantime20s_deadtime8min';
    %DataSet = 'Xenopus_atlas/Aug01_22-10_stage35p0_22p0keV_0700mm_50ms_1800proj';
    DataSet = 'Xenopus_atlas/Aug01_22-10_stage35p0_22p0keV_0700mm_50ms_1800proj';
end
if ~exist('FilterPixelThresholds','var')
    % Hot and dark pixel filter thresholds for dark fields, flat fields,
    % and sample images
    FilterPixelThresholds = [1.21 0.76,  1.035  0.97, 1.035 0.965];
end
if ~exist('PixelRegion','var')
    % PixelRegion = {ROWS, COLS}; ROWS = [START STOP]; COLS = [START STOP];
    %PixelRegion = [];%{[400 1200],[61 1400]};
    PixelRegion = {};%{[901 1028],[701 1212]};
end
if ~exist('EnergyDistancePixelsize','var')
    EnergyDistancePixelsize = [22.0 0.74 1.1e-6];
end
if ~exist('RegPar','var')
    RegPar = 2.5;
end
if ~exist('darkMeanThresh','var')
    darkMeanThresh = 1000;
end
if ~exist('PhaseMethod','var')
    PhaseMethod = 'quasi2';
end
if ~exist('BinaryFilterThreshold','var')
    BinaryFilterThreshold = 0.1;
end
if ~exist('saveIntMaps','var')
    % saves intensities as 'edf' or 'tif', or not at all
    saveIntMaps = '';'edf'; %'edf' or 'tif'
end
if ~exist('savePhaseMaps','var')
    savePhaseMaps = 'edf';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open pool of MATLAB sessions for parallel computation
if matlabpool('size') == 0
    matlabpool('open')
    mlpSize = matlabpool('SIZE');
else
    mlpSize = matlabpool('SIZE');
end
if mlpSize == 0
    mlpSize = 1;
end
%% Parameter
Energy    = EnergyDistancePixelsize(1);
Distance  = EnergyDistancePixelsize(2);
Pixelsize = EnergyDistancePixelsize(3);
lambda    = 6.62606896e-34*299792458/(Energy*1.60217733e-16);
ArgPrefac = 2*pi*lambda*Distance/Pixelsize^2;
hptDark = FilterPixelThresholds(1:2);
hptFlat = FilterPixelThresholds(3:4);
hptProj = FilterPixelThresholds(5:6);
projOffset = 1; % Dismiss first image of each acquisition sequence due to camera start up time
darkOffset = 4; % Dismiss first image of each acquisition sequence due to camera start up time
flatOffset = 4; % Dismiss first image of each acquisition sequence due to camera start up time
maxNumDarkFlat = 99; % Restrict maximum number of dark or flat fields to be used
%% Folder- and Filenames
CheckTrailingSlash(ParentPath);
CheckTrailingSlash(DataSet);
DataPath  = [ParentPath 'data/'];
InputPath = [DataPath DataSet];
IntPath   = sprintf('%s%s%s',ParentPath,'int/',DataSet);
PhasePath = sprintf('%s%s%s',ParentPath,'phase/',DataSet);
VolPath   = sprintf('%s%s%s',ParentPath,'vol/',DataSet);
%Prefix and postfix of output folder
Prefix = '';
Prefix = [Prefix 'FiltSino_'];
Prefix = [Prefix 'intMasking_'];
Postfix = '';
IntOutputFolder = 'int_filtSino';
fprintf('\nSTART DATA PROCESSING\nDATA SET: %s',DataSet(1:end-1))
% Number and indices of tomograms
TomoNames = FilenameCell(sprintf('%s*postDark_00001.tif',InputPath));
if isempty(TomoNames)
    fprintf('\n\n NO DARK FIELDS FOUND (Checked for ''postDark_00001.tif'')\n ABORT ''%s''\n\n',mfilename)
    return
end
NumTomosFound = numel(TomoNames);
%% Loop over tomograms
if isempty(TomoToProcess)
    for nn = NumTomosFound:-1:1
        a = TomoNames{nn}(6:8);
        b = regexp(a,'[^0-9]');
        a = a(1:(b(1)-1));
        a = str2double(a);
        TomoToProcess(nn) = a;
    end
    TomoToProcess = sort(TomoToProcess);
end
NumTomos = numel(TomoToProcess);
% Get image dimensions for preallocation of memory
if isempty(PixelRegion)
    filename = sprintf('%s%s',InputPath,TomoNames{1});
    iminfo = imfinfo(filename,'tif');
    dim1 = iminfo.Height;
    dim2 = iminfo.Width;
    PixelRegion = {[1 dim1],[1 dim2]};
else
    dim1 = PixelRegion{1}(2)-PixelRegion{1}(1)+1;
    dim2 = PixelRegion{2}(2)-PixelRegion{2}(1)+1;
end
fprintf('\n PROCESSING %u OF %u TOMOGRAMS: %s (index as it appears in the filename)',NumTomos,NumTomosFound,mat2str(TomoToProcess));
fprintf('\n HOT/DARK PIXEL FILTER THRESHOLDS: : [dark, flat, proj] = [%g %g, %g %g, %g %g]',FilterPixelThresholds);
fprintf('\n IMAGE DIMENSIONS: [ver x hor] = [%u x %u]',dim1,dim2);
fprintf('\n ROI: PixelRegion = {[%u %u],[%u %u]}',cell2mat(PixelRegion));
fprintf('\n REGULARIZATI0N: alpha = %g, 1/2/( 2*pi*lambda*z/dx^2 * xi^2/2 + 10^-alpha ) = 1/2/( %g * [-1/2,-1/2-1/N] + %g )',RegPar,ArgPrefac,10^-RegPar);
%% Phase retrieval filters
[~,PostfixBinFilt] = PhaseFilter(PhaseMethod,[dim1 dim2],EnergyDistancePixelsize,RegPar,BinaryFilterThreshold);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start loop over tomograms%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tuptoStartLoop = toc;
fprintf('\n Start looping over tomos. Elapsed time: %.2gs',tuptoStartLoop);
tloop = 0;
for tt = 1:NumTomos
    tic;
    tomoInd = TomoToProcess(tt);
    fprintf('\n\nPROCESSING TOMO NUMBER: %u, FILE INDEX: %u',tt,tomoInd)
    % Check number of flats, darks, and projections, and read in name patterns
    preDarkNames = FilenameCell(sprintf('%sproj_%upreDark_*.tif',InputPath,tomoInd));
    NumPreDark = numel(preDarkNames);
    preFlatNames = FilenameCell(sprintf('%sproj_%upreFlat_*.tif',InputPath,tomoInd));
    NumPreFlat = numel(preFlatNames);
    projNames = FilenameCell(sprintf('%sproj_%u_*.tif',InputPath,tomoInd));
    NumProj = numel(projNames);
    postFlatNames = FilenameCell(sprintf('%sproj_%upostFlat_*.tif',InputPath,tomoInd));
    NumPostFlat = numel(postFlatNames);
    postDarkNames = FilenameCell(sprintf('%sproj_%upostDark_*.tif',InputPath,tomoInd));
    NumPostDark = numel(postDarkNames);
    fprintf('\n IMAGES FOUND: %u preDark, %u preFlat, %u proj, %u postFlat, %u postDark', ...
        NumPreDark,NumPreFlat,NumProj,NumPostFlat,NumPostDark)
    if ~(((NumProj*NumPreDark*NumPreFlat>0) || (NumProj*NumPostFlat*NumPostFlat>0)))
        if tt < NumTomos
            fprintf('\n CONTINUE WITH NEXT LOOP INDEX! (TOMOGRAM NUMBER: %u. LOOP INDEX: %u.)\n',tomoInd,tt)
            continue;
        else
            fprintf('\n LOOP OVER TOMOGRAMS STOPPED! (TOMOGRAM NUMBER: %u. LOOP INDEX: %u.)\n',tomoInd,tt)
            return;
        end
    end
    % Check for superfluous projection which occured for some scans after
    % the last 'real' projection was acquired
    if NumProj < 800
        if NumProj > 500
            projNames(501:end) = [];
            NumProj = numel(projNames);
        end
    end
    dim3 = NumProj-projOffset;
    % Output paths and folder names
    if ~isempty(saveIntMaps)
        IntOutputPath = sprintf('%stomo%02u',IntPath,tomoInd);
        CheckTrailingSlash(IntOutputPath);
        CheckAndMakePath(IntOutputPath);
        fprintf('\n INT OUTPUT PATH: %s',IntOutputPath);
    end
    PhaseFolderName = sprintf('%sFDcor_%s_regPar%3.2f%s%s/',Prefix,PhaseMethod,RegPar,PostfixBinFilt,Postfix);
    PhaseFolderName = regexprep(PhaseFolderName,'\.','p');
    OutputFolder = sprintf('%stomo%02u/',PhaseFolderName,tomoInd);
    if ~isempty(savePhaseMaps)
        fprintf('\n PHASE OUTPUT FOLDER: %s',OutputFolder);
        PhaseOutputPath = sprintf('%s%s',PhasePath,OutputFolder);
        CheckAndMakePath(PhaseOutputPath);
    end
    VolOutputPath = sprintf('%s%s',VolPath,PhaseFolderName);
    CheckAndMakePath(VolOutputPath);
    %% PREDARK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    darkVariant = 0;
    if NumPreDark > 2
        NumPreDark = min(NumPreDark,maxNumDarkFlat);
        stack = zeros(dim1,dim2,NumPreDark-darkOffset,'single');
        parfor nn = 1:NumPreDark-darkOffset
            filename = sprintf('%s%s',InputPath,preDarkNames{nn+darkOffset});
            stack(:,:,nn) = FilterPixel(single(imread(filename,'tif','PixelRegion',PixelRegion)),hptDark);
        end
        preDarkMean = mean(stack,3);
        darkMeanAv = mean(preDarkMean(:));
        if darkMeanAv > darkMeanThresh
            fprintf('\n\n !!PREDARK FIELDS SEEM TO HAVE BEEN EXPOSED!! MEAN VALUE OF MEAN DARK: %g\n',darkMeanAv);
        else
            darkVariant = darkVariant + 1;
        end
    end
    %% POSTDARK
    if NumPostDark > 2
        NumPostDark = min(NumPostDark,maxNumDarkFlat);
        stack = zeros(dim1,dim2,NumPostDark-darkOffset,'single');
        parfor nn = 1:NumPostDark-darkOffset
            filename = sprintf('%s%s',InputPath,postDarkNames{nn+darkOffset});
            stack(:,:,nn) = FilterPixel(single(imread(filename,'tif','PixelRegion',PixelRegion)),hptDark);
        end
        postDarkMean = mean(stack,3);
        darkMeanAv = mean(postDarkMean(:));
        if darkMeanAv > darkMeanThresh
            fprintf('\n\n !!POSTDARK FIELDS SEEM TO HAVE BEEN EXPOSED!! MEAN VALUE OF MEAN DARK: %g\n',darkMeanAv);
        else
            darkVariant = darkVariant + 2;
        end
    end
    % Final dark field
    switch darkVariant
        case 0
            darkMean = 0;
            fprintf('\n DARK FIELD IS SET TO ZERO.')
        case 1
            darkMean = preDarkMean;
            fprintf('\n ONLY PREDARK FIELD IS USED.')
        case 2
            darkMean = postDarkMean;
            fprintf('\n ONLY POSTDARK FIELD IS USED.')
        case 3
            darkMean = (preDarkMean + postDarkMean)/2;
    end
    clear preDarkMean postDarkMean
    %% PREFLAT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    flatVariant = 0;
    if NumPreFlat > 2
        NumPreFlat = min(NumPreFlat,maxNumDarkFlat);
        stack = zeros(dim1,dim2,NumPreFlat-flatOffset,'single');
        parfor nn = 1:NumPreFlat-flatOffset
            filename = sprintf('%s%s',InputPath,preFlatNames{nn+flatOffset});
            stack(:,:,nn) = FilterPixel(single(imread(filename,'tif','PixelRegion',PixelRegion)),hptFlat);
        end
        preFlatMean = mean(stack,3);
        flatVariant = flatVariant + 1;
    end
    %% POSTFLAT
    if NumPostFlat > 2
        NumPostFlat = min(NumPostFlat,maxNumDarkFlat);
        stack = zeros(dim1,dim2,NumPostFlat-flatOffset,'single');
        parfor nn = 1:NumPostFlat-flatOffset
            filename = sprintf('%s%s',InputPath,postFlatNames{nn+flatOffset});
            stack(:,:,nn) = FilterPixel(single(imread(filename,'tif','PixelRegion',PixelRegion)),hptFlat);
        end
        postFlatMean = mean(stack,3);
        flatVariant = flatVariant + 2;
    end
    switch flatVariant
        case 0
            fprintf('\n !!NO FLAT FIELDS FOUND. EXITING LOOP!!\n\n')
            return
        case 1
            flatMean = preFlatMean;
            fprintf(' ONLY PREFLATS FOUND.')
        case 2
            flatMean = postFlatMean;
            fprintf(' ONLY POSTFLATS FOUND.')
        case 3
            flatMean = (preFlatMean + postFlatMean)/2;
            fprintf(' PRE- AND POSTFLATS FOUND.')
    end
    clear preFlatMean postFlatMean
    %% Modify dark field.
    % Needed because aperture of lens is seen on images
    darkMeanAv = mean(darkMean(:));
    mask = ones(dim1,dim2);
    mask(flatMean < 5*mean(darkMeanAv)) = 0;
    mask = medfilt2(mask,[3 3],'symmetric');
    %mask = imfilter(mask,fspecial('gaussian',[3 3],10),'symmetric');
    mask = imfilter(mask,fspecial('disk',10),'symmetric');
    darkMean = darkMean.*mask;
    flatMean = flatMean - darkMean;
    mask = floor(mask);
    %% ROTATION AXIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Read first projection, read last projection and flip
    filename = sprintf('%s%s',InputPath,projNames{1+projOffset});
    im = (FilterPixel(single(imread(filename,'tif','PixelRegion',PixelRegion)),hptProj) - darkMean)./flatMean;
    filename = sprintf('%s%s',InputPath,projNames{NumProj});
    im180 = (FilterPixel(single(imread(filename,'tif','PixelRegion',PixelRegion)),hptProj) - darkMean)./flatMean;
    % Cropping for determination of rot axis
    % Cropping values for determination of rotation axis
    xi = 0;
    xf = 0;
    % vertical
    leftCropIndForRotAx = ceil(dim2/4);
    if (PixelRegion{1}(1) < xi) && (PixelRegion{1}(2) > xf)
        x = (xi:xf)-PixelRegion{1}(1);
    else
        if dim1 > 300
            x = floor(3/10*dim1):ceil(dim1*(10-3)/10);
        else
            x = 1:dim1;
        end
    end
    % horizontal
    if (PixelRegion{2}(1) < leftCropIndForRotAx) && (PixelRegion{1}(2) > (dim2-leftCropIndForRotAx))
        leftCropIndForRotAx = leftCropIndForRotAx - PixelRegion{2}(1) + 1;
        y = leftCropIndForRotAx:dim2-(leftCropIndForRotAx-1);
    else
        if dim2 > 300
            leftCropIndForRotAx = ceil(3/10*dim2);
            y = leftCropIndForRotAx:dim2-(leftCropIndForRotAx-1);
        else
            y = 1:dim2;
            leftCropIndForRotAx = 1;
        end
    end
    im = im(x,y);
    im180 = im180(x,y);
    % Compute axis
    out        = ImageCorrelation(im,fliplr(im180),0,0);
    RotAxisPos = out.VerticalRotationAxisPosition + (leftCropIndForRotAx - 1);
    RotAxisPosPadded = RotAxisPos;
    fprintf('\n ROTATION AXIS: %6g',RotAxisPos);    
    clear im im180;
    fprintf('\n Rotation axis of ROI {[ver],[hor]} = {[%u %u],[%u %u]} of [ver hor] = [%u %u]',x(1),x(end),y(1),y(end),dim1,dim2);
    % Make par file
    if ~isempty(savePhaseMaps)
        ParInputFilePrefix = sprintf('%sphase_',PhaseOutputPath);
        ParOutputFilePrefix = sprintf('%stomo%02u.vol',VolOutputPath,tomoInd);
        MakeParFile32ID(ParInputFilePrefix,ParOutputFilePrefix,[dim1 dim2],RotAxisPosPadded,'NumberOfProjections',dim3,'OverallAngle',180*(1-projOffset/NumProj));
    end
%     if ~isempty(saveIntMaps)
%         ParInputFilePrefix = sprintf('%sint_',IntOutputPath);
%         IntVolOutputPath = sprintf('%s%s',VolPath,IntOutputFolder);
%         CheckTrailingSlash(IntVolOutputPath);
%         CheckAndMakePath(IntVolOutputPath);
%         ParOutputFilePrefix = sprintf('%stomo%02u.vol',IntVolOutputPath,tomoInd);
%         MakeParFile32ID(ParInputFilePrefix,ParOutputFilePrefix,[dim1 dim2],RotAxisPosPadded,'NumberOfProjections',dim3,'OverallAngle',180*(1-projOffset/NumProj));
%     end
    fprintf('\n Mean flat and mean dark computed. Elapsed time: %.2g min',toc/60);
    fprintf('\n Start: reading, pixel filtering, flat/dark correction');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PROJECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    stack = zeros(dim1,dim2,dim3,'single');
    %% Reading, hot-pixel filtering, flat/dark field, cropping/padding
    stackSize = Bytes(stack,0);
    if stackSize ~= dim1*dim2*dim3*4
        fprintf('\n\n ERROR: SIZE OF PROJECTION STACK (%u B) UNEQUAL DIM1*DIM2*DIM3 (%u B)!\n\n',stackSize,dim1*dim2*dim3*4)
        return;
    end
    % Limitation for parallel computing: The data which will be transfered
    % from server to client is limited to ~2 GB
    parFlag = stackSize / mlpSize < 2*1024^3;
    fprintf('\n Size of projection array: %.2g GB (%u B)',stackSize/1024^3,stackSize);
    fprintf('\n Parallel loop flag: %g',parFlag);
    if parFlag
        parfor nn = 1:dim3
            filename = sprintf('%s%s',InputPath,projNames{nn+projOffset});
            im = FilterPixel(single(imread(filename,'tif','PixelRegion',PixelRegion)),hptProj);
            im = (im - darkMean)./flatMean;
            im = im.*mask;
            im = im + sum(im(:))/sum(mask(:))*(1-mask);
            %             switch saveIntMaps
            %                 case 'edf'
            %                     edfwrite(sprintf('%sint_%04u.edf',IntOutputPath,nn),im','float32');
            %                 case 'tif'
            %                     write32bitTIFfromSingle(sprintf('%sint_%04u.tif',IntOutputPath,nn),im);
            %             end
            stack(:,:,nn) = im;
        end
    else
        fprintf('\n  READING PROJECTION NUMBER (stack size forbids usage of parallel computing):\n   ');
        for nn = dim3:-1:1
            PrintNum(nn);
            filename = sprintf('%s%s',InputPath,projNames{nn+projOffset});
            im = FilterPixel(single(imread(filename,'tif','PixelRegion',PixelRegion)),hptProj);
            im = (im - darkMean)./flatMean;
            im = im.*mask;
            im = im + sum(im(:))/sum(mask(:))*(1-mask);
            %             switch saveIntMaps
            %                 case 'edf'
            %                     edfwrite(sprintf('%sint_%04u.edf',IntOutputPath,nn),im','float32');
            %                 case 'tif'
            %                     write32bitTIFfromSingle(sprintf('%sint_%04u.tif',IntOutputPath,nn),im);
            %             end
            stack(:,:,nn) = im;
        end
    end
    fprintf('\n Elapsed time: %.2g min',toc/60);
    %% Phase retrieval and sinogram filtering
    %fprintf('\n Combined phase retrieval and sinogram filtering. Elapsed time: %.2g min',toc/60);
    %stack = FilterStack(stack,doFilterSino,PhaseMethod,EnergyDistancePixelsize,RegPar,BinaryFilterThreshold);
    %stack = PhaseSinoFilter(stack,PhaseMethod,EnergyDistancePixelsize,RegPar,BinaryFilterThreshold);
    % Fourier space filter for phase retrieval
    phaseFilter = PhaseFilter(PhaseMethod,size(stack),EnergyDistancePixelsize,RegPar,BinaryFilterThreshold);
    % Due Matlab's memory handling an memory limitation, sino filtering
    % and phase retrieval will not be combined
    %% Sinogram filtering
    fprintf('\n Sino filtering of sino slice:\n');
    for nn = 1:dim1
        PrintNum(nn);
        im = squeeze(stack(nn,:,:));
        im = fft2(im);
        im(:,1) = median(im(:,[1:3 end-1:end]),2);
        im = real(ifft2(im));
        stack(nn,:,:) = im;
    end
    fprintf('\n Elapsed time: %.2g min',toc/60);
    % Make .par file for intensities
    if ~isempty(saveIntMaps)
        ParInputFilePrefix = sprintf('%sint_',IntOutputPath);
        IntVolOutputPath = sprintf('%s%s',VolPath,IntOutputFolder);
        CheckTrailingSlash(IntVolOutputPath);
        CheckAndMakePath(IntVolOutputPath);
        ParOutputFilePrefix = sprintf('%stomo%02u.vol',IntVolOutputPath,tomoInd);
        MakeParFile32ID(ParInputFilePrefix,ParOutputFilePrefix,[dim1 dim2],RotAxisPosPadded,'NumberOfProjections',dim3,'OverallAngle',180*(1-projOffset/NumProj));
    end
    %% Phase retrieval
    fprintf('\n Phase retrieval of projection:\n');
    for nn = 1:dim3
        PrintNum(nn);
        im = squeeze(stack(:,:,nn));
        switch saveIntMaps
            case 'edf'
                edfwrite(sprintf('%sint_%04u.edf',IntOutputPath,nn),im','float32');
            case 'tif'
                write32bitTIFfromSingle(sprintf('%sint_%04u.tif',IntOutputPath,nn),im);
        end
        im = real(ifft2(phaseFilter.*fft2(im)));
        stack(:,:,nn) = im;
        switch savePhaseMaps
            case 'edf'
                edfwrite(sprintf('%sphase_%04u.edf',PhaseOutputPath,nn),im','float32');
            case 'tif'
                write32bitTIFfromSingle(sprintf('%sphase_%04u.tif',PhaseOutputPath,nn),im);
        end
    end
    fprintf('\n Sino filtering and phase retrieval finished.');
    fprintf('\n Elapsed time: %.2g min',toc/60);
    tloop = tloop + toc;
end
ttotal = tuptoStartLoop + tloop;
fprintf('\nTOTAL PROCESSING TIME: %g s = %.2g min = %.4g h\n\n',ttotal,ttotal/60,ttotal/60/60);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
end