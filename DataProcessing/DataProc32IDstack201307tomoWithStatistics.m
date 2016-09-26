function DataProc32IDstack201307tomo(varargin)
% data processing for the experiment In vivo imaging, GUP at beamline
% 32ID-C at APS, ANL, Chicago, Illinois, USA, 2013-07-28 to 2013-08-02.
% Allows for tomographic reconstruction. Hower, due to Matlab's slow
% interpolation function, the use of PyHST is much faster.
%
%Written by Julian Moosmann, last modified 2013-09-13

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
%% PARSE INPUT ARGUMENT 'VARARGIN'
ParseVarargin(varargin)
%% DEFAULT ARGUMENTS AND PARAMETERS
if ~exist('TomoToProcess','var')
    % Index of tomogram as it appears in the filename. It's not the
    % absolute number of tomograms
    TomoToProcess = [];
end
if ~exist('ParentPath','var')
    ParentPath = '/mnt/tomoraid-LSDF/tomo/APS_32ID-C_InVivoImaging_GUP34879_2013-07-30/';    
end
if ~exist('DataSet','var')
    %DataSet = 'Lamprey_atlas/Aug02_01-01_Lamprey_stage21p0_22p0keV_0700mm_50ms_1800proj';
    %DataSet = 'Lamprey_inVivo/August01_07-10_Lamprey_stage10somites__30p0keV_0700mm_15ms_500proj_scantime20s_deadtime08min';
    %DataSet = 'Xenopus_inVivo/Jul28_14-50_urea_stage17p0_30p0keV_0400mm_20ms_0500proj_scantime20s_deadtime8min';
    %DataSet = 'Zebrafish_inVivo/July31_20-00_Zebrafish_stage07p0hpf_30p0keV_0700mm_15ms_500proj_scantime20s_deadtime08min';
    %DataSet = 'Zebrafish_atlas/August02_02-30_Zebrafish_stage20p0hpf_22p0keV_0700mm_50ms_1800proj';    
    %DataSet = 'Xenopus_atlas/stage35p0';
    %DataSet = 'Xenopus_inVivo/Jul28_20-25_urea_stage21p0_30p0keV_0400mm_10ms_0500proj_scantime20s_deadtime8min';
    %DataSet = 'Xenopus_inVivo/Jul29_15-10_urea_stage27p0_30p0keV_0700mm_15ms_0500proj_scantime20s_deadtime8min';
    %DataSet = 'Xenopus_atlas/stage41p0';
    DataSet = 'Jul29_05-15_urea_stage25p0_30p0keV_0400mm_10ms_0500proj_scantime20s_deadtime8min';
end
if ~exist('writeFormat','var')
    writeFormat = 'edf';
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
    EnergyDistancePixelsize = [30.0 0.7 1.1e-6];
end
if ~exist('RegPar','var')
    RegPar = 2.5;
end
if ~exist('darkMeanThresh','var')
    darkMeanThresh = 1000;
end
if ~exist('PhaseMethod','var')
    PhaseMethod = 'tie';
end
if ~exist('BinaryFilterThreshold','var')
    BinaryFilterThreshold = 0.2;
end
if ~exist('doMeanSubstraction','var')
    doMeanSubstraction = 1;
end
if ~exist('doFilterSino','var')
    doFilterSino = 1;
end
if ~exist('doDarkFieldCorrection','var')
    doDarkFieldCorrection = 1;
end
if ~exist('doIntMasking','var')
    doIntMasking(1) = 1;
end
if ~exist('saveIntMaps','var')
    % saves intensities as 'edf' or 'tif', or not at all 
    saveIntMaps = 'edf'; %'edf' or 'tif'
end

if ~exist('savePhaseMaps','var')
    savePhaseMaps(1) = 0;
end
if ~exist('doFBP','var')
    doFBP(1) = 0;
end
if doFBP
    if ~exist('iradonInterpolation','var')
        iradonInterpolation = 'spline';% 'nearest' 'spline' 'pchip' 'cubic' 'v5cubic'
    end
    if ~exist('iradonFilter','var')
        iradonFilter = 'Ram-Lak';% 'Shepp-Logan' 'Cosine' 'Hamming' 'Hann'  'none'
    end
    if ~exist('iradonFrequencyScaling','var')
        iradonFrequencyScaling = 1; %  scalar in (0,1]. If FREQUENCY_SCALING is less than 1, the filter is compressed to fit into the frequency range [0,FREQUENCY_SCALING], in normalized frequencies; all frequencies above FREQUENCY_SCALING are set to 0.
    end
    if ~exist('iradonOutputSize','var')
        iradonOutputSize = 1448;
        %iradonOutputSize = 2*floor(size(R,1)/(2*sqrt(2))); this square
        %fits into the circle for which there is complete information
        %available from the projections
    end
end
if ~exist('showFigures','var')
    showFigures(1) = 0;
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
projOffset = 2; % Dismiss first image of each acquisition sequence due to camera start up time
darkOffset = 4; % Dismiss first image of each acquisition sequence due to camera start up time
flatOffset = 4; % Dismiss first image of each acquisition sequence due to camera start up time
maxNumDarkFlat = 100; % Restrict maximum number of dark or flat fields to be used
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
Prefix = [Prefix '3DstackProc_'];
if ~doDarkFieldCorrection,Prefix = [Prefix 'noDark_'];end
if doFilterSino,Prefix = [Prefix 'FiltSino_'];end
if doIntMasking,Prefix = [Prefix 'intMasking_'];end
Postfix = '';
if ~doMeanSubstraction,Postfix = '_noMeanSub';end
IntOutputFolder = 'int';
% if doIntMasking
%     IntOutputFolder = [IntOutpuFolder '_withMasking'];
% end    
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
% Preallocation
hpfPreDarkAr = cell(NumTomos,1);
hpfPostDarkAr = cell(NumTomos,1);
hpfPreFlatAr = cell(NumTomos,1);
hpfPostFlatAr = cell(NumTomos,1);
hpfProjAr = cell(NumTomos,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start loop over tomograms%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tuptoStartLoop = toc;
fprintf('\n Time up to loop over tomo: %.2gs',tuptoStartLoop);
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
    if savePhaseMaps
        fprintf('\n PHASE OUTPUT FOLDER: %s',OutputFolder);
        PhaseOutputPath = sprintf('%s%s',PhasePath,OutputFolder);
        CheckAndMakePath(PhaseOutputPath);
    end    
    VolOutputPath = sprintf('%s%s',VolPath,PhaseFolderName);
    CheckAndMakePath(VolOutputPath);
    if doFBP
        SlicesOutputPath = sprintf('%stomo%02u_%s_%s',VolOutputPath,tomoInd,iradonInterpolation,writeFormat(1:3));
        CheckTrailingSlash(SlicesOutputPath);
        CheckAndMakePath(SlicesOutputPath);
    end
    %% PREDARK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    darkVariant = 0;
    if NumPreDark > 2
        NumPreDark = mod(NumPreDark,maxNumDarkFlat);
        stack = zeros(dim1,dim2,NumPreDark-darkOffset,'single');
        hpfPreDark = zeros(1,size(stack,3),'single');
        parfor nn = 1:NumPreDark-darkOffset
            filename = sprintf('%s%s',InputPath,preDarkNames{nn+darkOffset});
            [stack(:,:,nn), hpfPreDark(nn)] = FilterPixel(single(imread(filename,'tif','PixelRegion',PixelRegion)),hptDark);
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
        NumPostDark = mod(NumPostDark,maxNumDarkFlat);
        stack = zeros(dim1,dim2,NumPostDark-darkOffset,'single');
        hpfPostDark = zeros(1,size(stack,3),'single');
        parfor nn = 1:NumPostDark-darkOffset
            filename = sprintf('%s%s',InputPath,postDarkNames{nn+darkOffset});
            [stack(:,:,nn), hpfPostDark(nn)] = FilterPixel(single(imread(filename,'tif','PixelRegion',PixelRegion)),hptDark);
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
    if ~doDarkFieldCorrection
        darkVariant = 0;
    end
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
        NumPreFlat = mod(NumPreFlat,maxNumDarkFlat);
        stack = zeros(dim1,dim2,NumPreFlat-flatOffset,'single');
        hpfPreFlat = zeros(1,size(stack,3),'single');
        parfor nn = 1:NumPreFlat-flatOffset
            filename = sprintf('%s%s',InputPath,preFlatNames{nn+flatOffset});
            [stack(:,:,nn), hpfPreFlat(nn)] = FilterPixel(single(imread(filename,'tif','PixelRegion',PixelRegion)),hptFlat);
        end
        preFlatMean = mean(stack,3);
        flatVariant = flatVariant + 1;
    end
    %% POSTFLAT
    if NumPostFlat > 2
        NumPostFlat = mod(NumPostFlat,maxNumDarkFlat);
        stack = zeros(dim1,dim2,NumPostFlat-flatOffset,'single');
        hpfPostFlat = zeros(1,size(stack,3),'single');
        parfor nn = 1:NumPostFlat-flatOffset
            filename = sprintf('%s%s',InputPath,postFlatNames{nn+flatOffset});
            [stack(:,:,nn), hpfPostFlat(nn)] = FilterPixel(single(imread(filename,'tif','PixelRegion',PixelRegion)),hptFlat);
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
    if doDarkFieldCorrection
        darkMeanAv = mean(darkMean(:));
        mask = ones(dim1,dim2);
        mask(flatMean < 5*mean(darkMeanAv)) = 0;
        mask = medfilt2(mask,[3 3],'symmetric');
        %mask = imfilter(mask,fspecial('gaussian',[3 3],10),'symmetric');
        mask = imfilter(mask,fspecial('disk',10),'symmetric');
        darkMean = darkMean.*mask;
        flatMean = flatMean - darkMean;
        if doIntMasking
            mask = floor(mask);
        else
            clear mask;
        end
    end
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
    if doFBP
        % Crop/pad images to a power of who such that center coincides with rotation axis as iradon assumes
        im  = RotAxisSymmetricPadding(im,RotAxisPos);
        im180 = RotAxisSymmetricPadding(im180,RotAxisPos);
        out = ImageCorrelation(im,fliplr(im180),0,0);
        RotAxisPosPadded = out.VerticalRotationAxisPosition;
        % !!!! REASSIGN A NEW VALUE TO HORIZONTAL DIMENSION VARIABLE 'dim2' !!!
        newdim2 = size(im,2);
        fprintf('\n NEW IMAGE DIMENSIONS: [ver x hor] = [%u x %u]',dim1,newdim2)
        fprintf('\n ROTATION AXIS BEFORE PAD/CROP: %6g',RotAxisPos);
        fprintf('\n ROTATION AXIS AFTER  PAD/CROP: %6g',RotAxisPosPadded);
    else
        newdim2 = dim2;
        RotAxisPosPadded = RotAxisPos;
         fprintf('\n ROTATION AXIS: %6g',RotAxisPos);
    end
    clear im im180;
    fprintf('\n Rotation axis of ROI {[ver],[hor]} = {[%u %u],[%u %u]} of [ver hor] = [%u %u]',x(1),x(end),y(1),y(end),dim1,dim2);
    % Make par file
    if savePhaseMaps
        ParInputFilePrefix = sprintf('%sphase_',PhaseOutputPath);
        ParOutputFilePrefix = sprintf('%stomo%02u.vol',VolOutputPath,tomoInd);
        MakeParFile32ID(ParInputFilePrefix,ParOutputFilePrefix,[dim1 newdim2],RotAxisPosPadded,'NumberOfProjections',dim3,'OverallAngle',180*(1-projOffset/NumProj));
    end
    if ~isempty(saveIntMaps)
        ParInputFilePrefix = sprintf('%sint_',IntOutputPath);
        IntVolOutputPath = sprintf('%s%s',VolPath,IntOutputFolder);
        CheckTrailingSlash(IntVolOutputPath);
        CheckAndMakePath(IntVolOutputPath);
        ParOutputFilePrefix = sprintf('%stomo%02u.vol',IntVolOutputPath,tomoInd);
        MakeParFile32ID(ParInputFilePrefix,ParOutputFilePrefix,[dim1 newdim2],RotAxisPosPadded,'NumberOfProjections',dim3,'OverallAngle',180*(1-projOffset/NumProj));
    end
    fprintf('\n Mean flat and mean dark computed. Elapsed time: %.2g min',toc/60);
    fprintf('\n Start: reading, pixel filtering, flat/dark correction');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PROJECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if doFBP || savePhaseMaps
        stack = zeros(dim1,newdim2,dim3,'single');
        hpfProj = zeros(1,dim3,'single');
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
                [im,hpfProj(nn)] = FilterPixel(single(imread(filename,'tif','PixelRegion',PixelRegion)),hptProj);
                im = (im - darkMean)./flatMean;
                if doIntMasking && doDarkFieldCorrection
                    im = im.*mask;
                    im = im + sum(im(:))/sum(mask(:))*(1-mask);
                end
                im = RotAxisSymmetricPadding(im,RotAxisPos*doFBP);
                switch saveIntMaps
                    case 'edf'
                        edfwrite(sprintf('%sint_%04u.edf',IntOutputPath,nn),im','float32');
                    case 'tif'
                        write32bitTIFfromSingle(sprintf('%sint_%04u.tif',IntOutputPath,nn),im);
                end
                stack(:,:,nn) = im;
            end
        else
            fprintf('\n  READING PROJECTION NUMBER (stack size forbids usage of parallel computing):\n   ');
            for nn = 1:dim3
                fprintf('%4u ',nn);
                if mod(nn,20)==0
                    fprintf('\n   ');
                end
                filename = sprintf('%s%s',InputPath,projNames{nn+projOffset});
                [im,hpfProj(nn)] = FilterPixel(single(imread(filename,'tif','PixelRegion',PixelRegion)),hptProj);
                im = (im - darkMean)./flatMean;
                if doIntMasking
                    im = im.*mask;
                    im = im + sum(im(:))/sum(mask(:))*(1-mask);
                end
                im = RotAxisSymmetricPadding(im,RotAxisPos*doFBP);
                switch saveIntMaps
                    case 'edf'
                        edfwrite(sprintf('%sint_%04u.edf',IntOutputPath,nn),im','float32');
                    case 'tif'
                        write32bitTIFfromSingle(sprintf('%sint_%04u.tif',IntOutputPath,nn),im);
                end
                stack(:,:,nn) = im;
            end
        end
    else
        parfor nn = 1:dim3
            filename = sprintf('%s%s',InputPath,projNames{nn+projOffset});
            [im,hpfProj(nn)] = FilterPixel(single(imread(filename,'tif','PixelRegion',PixelRegion)),hptProj);
            im = (im - darkMean)./flatMean;
            if doIntMasking && doDarkFieldCorrection
                im = im.*mask;
                im = im + sum(im(:))/sum(mask(:))*(1-mask);
            end
            switch saveIntMaps
                case 'edf'
                    edfwrite(sprintf('%sint_%04u.edf',IntOutputPath,nn),im','float32');
                case 'tif'
                    write32bitTIFfromSingle(sprintf('%sint_%04u.tif',IntOutputPath,nn),im);
            end            
        end
    end
    %% Combined phase retrieval and sinogram filtering
    if doFBP || savePhaseMaps
        stack = FilterStack(stack,doFilterSino,PhaseMethod,EnergyDistancePixelsize,RegPar,BinaryFilterThreshold);
    end
    %% Save phase maps
    if savePhaseMaps
        switch parFlag
            case 0
                switch writeFormat
                    case 'edf'
                        fprintf('\n SAVING PROJECTION NUMBER:\n  ');
                        for nn = 1:dim3
                            fprintf('%4u ',nn);
                            if mod(nn,20)==0
                                fprintf('\n');
                            end
                            %saved phases have correct orientation. edf header: Dim_1=horizontal, Dim_2=vertical
                            edfwrite(sprintf('%sphase_%04u.edf',PhaseOutputPath,nn),squeeze(stack(:,:,nn))','float32');
                        end
                    case {'tif','tiff'}
                        for nn= 1:dim3
                            write32bitTIFfromSingle(sprintf('%sphase_%04u.tif',PhaseOutputPath,nn),squeeze(stack(:,:,nn)));
                        end
                end
            case 1
                switch writeFormat
                    case 'edf'
                        parfor nn = 1:dim3
                            %saved phases have correct orientation. edf header: Dim_1=horizontal, Dim_2=vertical
                            edfwrite(sprintf('%sphase_%04u.edf',PhaseOutputPath,nn),squeeze(stack(:,:,nn))','float32');
                            % save sinograms !! INDEX NEEDS TO BE CHANGED
                            %edfwrite(sprintf('%ssino_phase_z%04u.edf',PhaseOutputPath,nn),squeeze(stack(nn,:,:)),'float32');
                        end
                    case {'tif','tiff'}
                        parfor nn= 1:dim3
                            write32bitTIFfromSingle(sprintf('%sphase_%04u.tif',PhaseOutputPath,nn),squeeze(stack(:,:,nn)));
                        end
                end
        end
    end
    %% Inverse radon transformation
    if doFBP
        fprintf('\n Start iradon for tomographic reconstruction. Elapsed time: %.2g min',toc/60);
        theta = 180/NumProj; % theta = angles in degrees or angular increment
        parfor nn = 1:dim1 % in Matlab notation dim1 corresponds to the vertical direction
            im = iradon2(squeeze(stack(nn,:,:)),theta,iradonInterpolation,iradonFilter,iradonFrequencyScaling,iradonOutputSize);
            switch writeFormat
                case 'edf'
                    edfwrite(sprintf('%sphase_z%04u.edf',SlicesOutputPath,nn),im','float32');
                case {'tif','tiff'}
                    write32bitTIFfromSingle(sprintf('%sphase_z%04u.tif',SlicesOutputPath,nn),im);
            end
        end
        fprintf('\n iradon finished. Elapsed time: %.2 gmin',toc/60);
    else
        fprintf('\n Elapsed time: %.2g min',toc/60);
    end
    Numhpf = 0;
    if showFigures
        %% Gather hot pixel statistics
        
        if exist('hpfPreDark','var')
            hpfPreDarkAr{tt} = hpfPreDark;
            Numhpf = Numhpf + 1;
        end
        if exist('hpfPostDark','var')
            hpfPostDarkAr{tt} = hpfPostDark;
            Numhpf = Numhpf + 1;
        end
        if exist('hpfPreFlat','var')
            hpfPreFlatAr{tt} = hpfPreFlat;
            Numhpf = Numhpf + 1;
        end
        if exist('hpfPostFlat','var')
            hpfPostFlatAr{tt} = hpfPostFlat;
            Numhpf = Numhpf + 1;
        end
        if exist('hpfProj','var')
            hpfProjAr{tt}    = hpfProj;
            Numhpf = Numhpf + 1;
        end
    end
    tloop = tloop + toc;
end
ttotal = tuptoStartLoop + tloop;
fprintf('\nTOTAL PROCESSING TIME: %g s = %.2g min = %.4g h\n\n',ttotal,ttotal/60,ttotal/60/60);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Show hot pixel statistics
if showFigures*Numhpf > 0
    figure('Name','Percentage of filtered hot pixels')
    hpfInd = 1;
    if exist('hpfPreDark','var')
        subplot(Numhpf,1,hpfInd), plot(100*[hpfPreDarkAr{:}],'x'), title('PreDark')
        hpfInd = hpfInd + 1;
    end
    if exist('hpfPostDark','var')
        subplot(Numhpf,1,hpfInd), plot(100*[hpfPostDarkAr{:}],'x'), title('PostDark')
        hpfInd = hpfInd + 1;
    end
    if exist('hpfPreFlat','var')
        subplot(Numhpf,1,hpfInd), plot(100*[hpfPreFlatAr{:}],'x'), title('PreFlat')
        hpfInd = hpfInd + 1;
    end
    if exist('hpfPostFlat','var')
        subplot(Numhpf,1,hpfInd), plot(100*[hpfPostFlatAr{:}],'x'), title('PostFlat')
        hpfInd = hpfInd + 1;
    end
    if exist('hpfProj','var')
        subplot(Numhpf,1,hpfInd), plot(100*[hpfProjAr{:}],'x'), title('Projections')
    end
end
if Numhpf > 0
    save([PhasePath 'hpf.mat'],'hpf*Ar');
end
clear all
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stack = FilterStack(stack,doFilterSino,PhaseMethod,EnergyDistancePixelsize,RegPar,BinaryFilterThreshold)
% Child function for combined phase retrieval and sinogram filtering
% adapted for large array when memory limitations have to be taken into
% account..
maxSizeArrayGB = 40;
[dim1,dim2,dim3] = size(stack);
% Phase retrieval and sinogram filtering of 3D sinogram
fprintf('\n Compute FFT [2]. Elapsed time: %.2fmin. Stack size %.2g GB',toc/60,Bytes(stack,3))
stack = fft(stack,[],2);
%% Sino filtering
if doFilterSino
    fprintf('\n Start sino filtering. Elapsed time: %.2fmin. Stack size %.2g GB',toc/60,Bytes(stack,3))
    if Bytes(stack,3) < maxSizeArrayGB
        fprintf('\n Use direct 3D array manipulation: Apply FFT on whole stack');
        stack = fft(stack,[],3);
        fprintf('\n FFT for sino filtering finished. Elapsed time: %.2fmin. Stack size %.2g GB',toc/60,Bytes(stack,3))
        stack(:,:,1) = median(stack(:,:,[1:3 end-1:end]),3);
        fprintf('\n Median computation of slab finished. Elapsed time: %.2fmin',toc/60)
        stack = ifft(stack,[],3);
        fprintf('\n iFFT for sino filtering finished. Elapsed time: %.2fmin',toc/60)
    else
        fprintf('\n Use ''for'' loop instead of direct application of FT on 3D array because of array size: %.2g GB',Bytes(stack,3));
        fprintf('\n Sino filtering of projection number:\n  ');
        for nn = 1:size(stack,1)
            fprintf('%4u ',nn);if mod(nn,20)==0,fprintf('  \n   ');end
            sino = stack(nn,:,:);
            sino = fft(sino,[],3);
            sino(:,:,1) = median(sino(:,:,[1:3 end-1:end]),3);
            sino = ifft(sino,[],3);
            stack(nn,:,:) = sino;
        end
        fprintf('\n Sino filtering finished. Elapsed time: %.2fmin. Stack size %.2g GB',toc/60,Bytes(stack,3))
    end
end
stackSize = Bytes(stack,3);
fprintf('\n Compute FFT [1].')
if stackSize < maxSizeArrayGB
    stack = fft(stack,[],1);
else
    fprintf('\n Use of ''for'' loop over 3rd index n:\n  ');
    for nn = 1:dim3
        fprintf('%4u ',nn);if mod(nn,20)==0, fprintf('  \n   ');end
        stack(:,:,nn) = fft(stack(:,:,nn),[],1);
    end  
end
%% Phase retrieval
phaseFilter = PhaseFilter(PhaseMethod,[dim1 dim2],EnergyDistancePixelsize,RegPar,BinaryFilterThreshold);
stackSize = Bytes(stack,3);
fprintf('\n Apply phase filter. Elapsed time: %.2fmin. Stack size %.2g GB',toc/60,stackSize)
if stackSize < 30;
    stack = stack.*repmat(phaseFilter,[1 1 dim3]);
else
    fprintf('\n Use ''bsxfun'' for phase filter multiplication because of better memory efficency w.r.t. ''repmat''.')
    stack = bsxfun(@times,stack,phaseFilter);
end
fprintf('\n Compute inverse FFT [1,2]. Elapsed time: %.2 fmin',toc/60)
if stackSize < maxSizeArrayGB
    stack = ifft2(stack);%ifft(ifft(stack,[],1),[],2);
else
    fprintf('\n Use of ''for'' loop over 3rd index n:\n  ');
    for nn = 1:dim3
        fprintf('%4u ',nn);if mod(nn,20)==0,fprintf('  \n   ');end
        stack(:,:,nn) = ifft2(stack(:,:,nn));
    end
end
stack = real(stack);
fprintf('\n Phase retrieval and optional sino filtering finished. Elapsed time: %.2fmin',toc/60)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% END OF PRORGRAMME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
