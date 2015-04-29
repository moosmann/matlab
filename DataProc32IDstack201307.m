function DataProc32IDstack201307(TomoSetsToProcess,useGPU,doDarkFieldCorrection,LineFiltering,doFilterSino,doMeanSubstraction,doNormalizeSlices,showFigures,varargin)
% Data preprocessing for the experiment In vivo imaging, GUP at 32ID-C@APS

tic; tpretomoloop = 0; tpreprocloop = 0; tproc = 0; ttotal = 0;
%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    % Index of tomogram as it appears in the filename. It's not the
    % absolute number of tomograms
    TomoSetsToProcess = 1;%[];%
end
if nargin < 2
    useGPU = 0;
end
if nargin < 3
    doDarkFieldCorrection = 1;
end
if nargin < 4
    %[PreFDcorHor PostFDcorHor PostFDcorVert]
    LineFiltering = [0 0 0];
end
if nargin < 5
    doFilterSino = 0;
end
if nargin < 6
    doMeanSubstraction = 1;
end
if nargin < 7
    doNormalizeSlices = 0;
end
if nargin < 8
    showFigures = 1;
end
%% Input arguments
for nn = 1:2:numel(varargin)
    evalc([varargin{nn} ' = ' mat2str(varargin{nn+1})]);
    %assignin('caller',varargin{nn},varargin{nn+1});
end
%% Set default parameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('ParentPath','var')
    ParentPath = '/media/My Passport/APS_32ID-C_InVivoImaging_GUP34879_2013-07-30';
    
end
if ~exist('DataSet','var')
    %DataSet = 'Lamprey_atlas/Aug02_01-01_Lamprey_stage21p0_22p0keV_0700mm_50ms_1800proj';
    DataSet = 'Lamprey_inVivo/August01_07-10_Lamprey_stage10somites__30p0keV_0700mm_15ms_500proj_scantime20s_deadtime08min';
end
if ~exist('HotPixelThreshold','var')
    HotPixelThreshold = [1.2 1.1 1.1];
end
if ~exist('FSfiltering','var')
    FSfiltering = [0 0 0];
end
if ~exist('PixelRegion','var')
    % PixelRegion = {ROWS, COLS};
    % ROWS = [START STOP]; COLS = [START STOP];
    %PixelRegion = {};%{[400 1200],[61 1400]};
    PixelRegion = {[959 1058],[301 2016-300]};
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
    %PhaseMethod = 'pctf';
    %PhaseMethod = 'pctfHalfSine';
end
if ~exist('BinaryFilterThreshold','var')
    BinaryFilterThreshold = 0.2;
end
if ~exist('PipeGPUtag','var')
    PipeGPUtag = 0;
end
% if ~exist('doPreFilteringHorizontal','var')
%     doPreFilteringHorizontal = 0;
% end
% if ~exist('doPostFilteringHorizontal','var')
%     doPostFilteringHorizontal = 1;
% end
% if ~exist('TomoSetsToProcess','var')
%     TomoSetsToProcess = 3;%:13;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Open pool of MATLAB sessions for parallel computation
if useGPU
    if matlabpool('size') == 0
        matlabpool('open')
    end
end
%% Parameter
Energy    = EnergyDistancePixelsize(1);
Distance  = EnergyDistancePixelsize(2);
Pixelsize = EnergyDistancePixelsize(3);
lambda    = 6.62606896e-34*299792458/(Energy*1.60217733e-16);
ArgPrefac = 2*pi*lambda*Distance/Pixelsize^2;
hptDark = HotPixelThreshold(1);
hptFlat = HotPixelThreshold(2);
hptProj = HotPixelThreshold(3);
offset = 2; % Dismiss first image of each acquisition sequence due to camera start up time
maxNumDarkFlat = 100; % Restrict maximum number of dark or flat fields to be used
%% Paths
%PrePath(ParentPath)
CheckTrailingSlash(ParentPath);
CheckTrailingSlash(DataSet);
DataPath = [ParentPath 'data/'];
InputPath = [DataPath DataSet];
%IntPath  = [ParentPath 'int/'];
PhasePath  = [ParentPath 'phase/'];
CheckAndMakePath(PhasePath);
VolPath  = [ParentPath 'vol/'];
CheckAndMakePath(VolPath);
ParentOutputPath = sprintf('%s%s',PhasePath,DataSet);
ParentVolPath = sprintf('%s%s',VolPath,DataSet);
%Prefix and postfix of output folder
doPreFilteringHorizontal = LineFiltering(1);
doPostFilteringHorizontal = LineFiltering(2);
doPostFilteringVertical = LineFiltering(3);
Prefix = '';
if ~useGPU
    Prefix = [Prefix '3DstackProc_'];
end
if ~doDarkFieldCorrection
   Prefix = [Prefix 'noDark_']; 
end
if doPreFilteringHorizontal
   Prefix = [Prefix 'preFiltLineHor_']; 
end
if doPostFilteringHorizontal
    Prefix = [Prefix 'postFiltLineHor_'];
end
if doPostFilteringVertical
    Prefix = [Prefix 'postFiltLineVert_']; 
end
if doFilterSino && ~useGPU
    Prefix = [Prefix 'FiltSino_'];
end
Postfix = '';
if ~doMeanSubstraction
    Postfix = '_noMeanSub';
end
if doNormalizeSlices
    Postfix = [Postfix '_normSlices'];  
end
fprintf('\nSTARTING DATA PROCESSING\nDATA SET: %s\n',DataSet(1:end-1))
%fprintf(' DATA PATH: %s\n',DataPath)
% Number of tomograms found
TomoStruct = dir(sprintf('%s*postDark_00001.tif',InputPath));
if isempty(TomoStruct)
    fprintf('\n\n !!NO DARK FIELDS FOUND!! (Checked for ''postDark_00001.tif'')\n !!ABORT ''DataProc32ID''!!\n\n')
    return
end
NumTomosFound = numel(TomoStruct);
%% Loop over tomograms
if isempty(TomoSetsToProcess)
    for nn = 1:NumTomosFound
        a = TomoStruct(nn).name(6:8);
        b = regexp(a,'[^0-9]');
        a = a(1:(b(1)-1));
        a = str2double(a);
        TomoSetsToProcess(nn) = a;
    end
    TomoSetsToProcess = sort(TomoSetsToProcess);
end
NumTomos = numel(TomoSetsToProcess);
% Get image dimensions for preallocation of memory
if isempty(PixelRegion)
    filename = sprintf('%s%s',InputPath,TomoStruct(1).name);
    iminfo = imfinfo(filename,'tif');
    dim1 = iminfo.Height;
    dim2 = iminfo.Width;
    PixelRegion = {[1 dim1],[1 dim2]};
else
    dim1 = PixelRegion{1}(2)-PixelRegion{1}(1)+1;
    dim2 = PixelRegion{2}(2)-PixelRegion{2}(1)+1;
end
fprintf(' PROCESSING %u OF %u TOMOGRAMS: %s\n',NumTomos,NumTomosFound,mat2str(TomoSetsToProcess))
fprintf(' HOT PIXEL FILTER THRESHOLOD: : [dark flat proj] = [%g %g %g]\n',HotPixelThreshold)
fprintf(' IMAGE DIMENSIONS: [ver x hor] = [%u x %u]\n',dim1,dim2)
fprintf(' REGULARIZATI0N PARAMETER: alpha = %g, 1/2/( 2*pi*lambda*z/dx^2 * xi^2/2 + 10^-alpha ) = 1/2/( %g * [-1/2,-1/2-1/N] + %g )\n',RegPar,ArgPrefac,10^-RegPar)
%% Phase retrieval filters
[phaseFilter PostfixBinFilt] = PhaseFilter(PhaseMethod,[dim1 dim2],EnergyDistancePixelsize,RegPar,BinaryFilterThreshold);
% Preallocation
hpfPreDarkAr = cell(NumTomos,1);
hpfPostDarkAr = cell(NumTomos,1);
hpfPreFlatAr = cell(NumTomos,1);
hpfPostFlatAr = cell(NumTomos,1);
hpfProjAr    = cell(NumTomos,1);
tpretomoloop = tpretomoloop + toc;tic
%% Start looping %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for tt = 1:NumTomos
    tproc = 0;
    tomoInd = TomoSetsToProcess(tt);
    fprintf('\nPROCESSING TOMO INDEX: %u\n',tomoInd)
    %% Check number of flats, darks, and projections, and read in name patterns
    preDarkStruct = dir(sprintf('%sproj_%upreDark_*.tif',InputPath,tomoInd));
    NumPreDark = numel(preDarkStruct);
    preFlatStruct = dir(sprintf('%sproj_%upreFlat_*.tif',InputPath,tomoInd));
    NumPreFlat = numel(preFlatStruct);
    projStruct = dir(sprintf('%sproj_%u_*.tif',InputPath,tomoInd));
    NumProj = numel(projStruct);
    postFlatStruct = dir(sprintf('%sproj_%upostFlat_*.tif',InputPath,tomoInd));
    NumPostFlat = numel(postFlatStruct);
    postDarkStruct = dir(sprintf('%sproj_%upostDark_*.tif',InputPath,tomoInd));
    NumPostDark = numel(postDarkStruct);
    fprintf(' IMAGES FOUND: %u preDark, %u preFlat, %u proj, %u postFlat, %u postDark\n', ...
        NumPreDark,NumPreFlat,NumProj,NumPostFlat,NumPostDark)
    if ~(((NumProj*NumPreDark*NumPreFlat>0) || (NumProj*NumPostFlat*NumPostFlat>0)))
        fprintf('\n\n !!LOOP OVER TOMOGRAMS STOPPED! (TOMOGRAM NUMBER: %u. LOOP INDEX: %u.)\n\n',tomoInd,tt)
        return
    end
    % Check for superfluous projection which occured for some scans after
    % the actual tomogram
    if NumProj < 800
        NumProj = 500;
    end
    dim3 = NumProj-offset;
    %% Output paths and folder names
    PhaseFolderName = sprintf('%sFDcor_%s_regPar%3.2f%s%s/',Prefix,PhaseMethod,RegPar,PostfixBinFilt,Postfix);
    PhaseFolderName = regexprep(PhaseFolderName,'\.','p');
    OutputFolder = sprintf('%stomo%02u/',PhaseFolderName,tomoInd);
    fprintf(' OUTPUT FOLDER: %s\n',OutputFolder);
    OutputPath = sprintf('%s%s',ParentOutputPath,OutputFolder);
    CheckAndMakePath(OutputPath);
    VolOutputPath = sprintf('%s%s',ParentVolPath,PhaseFolderName);
    CheckAndMakePath(VolOutputPath);
    tpreprocloop = tpreprocloop + toc; tic
    fprintf(' PROCESSING: ');
    %% PREDARK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    darkVariant = 0;
    if NumPreDark > 2
        NumPreDark = mod(NumPreDark,maxNumDarkFlat);
        stack = zeros(dim1,dim2,NumPreDark-offset,'single');
        hpfPreDark = zeros(1,size(stack,3),'single');
        parfor nn = 1:NumPreDark-offset
            filename = sprintf('%s%s',InputPath,preDarkStruct(nn+offset).name);
            [stack(:,:,nn) hpfPreDark(nn)] = FilterHotPixel(single(imread(filename,'tif','PixelRegion',PixelRegion)),hptDark);
        end
        preDarkMean = mean(stack,3);
        darkMeanAv = mean(preDarkMean(:));
        if darkMeanAv > darkMeanThresh
            fprintf('\n\n !!PREDARK FIELDS SEEM TO HAVE BEEN EXPOSED!! MEAN VALUE OF MEAN DARK: %g\n\n PROCESSING: ',darkMeanAv);
        else
            darkVariant = darkVariant + 1;
        end
        t = toc; tproc = tproc + t; tic;
        fprintf('PreDark %gs. ',t)
    end
    %% POSTDARK
    if NumPostDark > 2
        NumPostDark = mod(NumPostDark,maxNumDarkFlat);
        stack = zeros(dim1,dim2,NumPostDark-offset,'single');
        hpfPostDark = zeros(1,size(stack,3),'single');
        parfor nn = 1:NumPostDark-offset
            filename = sprintf('%s%s',InputPath,postDarkStruct(nn+offset).name);
            [stack(:,:,nn) hpfPostDark(nn)] = FilterHotPixel(single(imread(filename,'tif','PixelRegion',PixelRegion)),hptDark);
        end
        postDarkMean = mean(stack,3);
        darkMeanAv = mean(postDarkMean(:));
        if darkMeanAv > darkMeanThresh
            fprintf('\n\n !!POSTDARK FIELDS SEEM TO HAVE BEEN EXPOSED!! MEAN VALUE OF MEAN DARK: %g\n\n PROCESSING: ',darkMeanAv);
        else
            darkVariant = darkVariant + 2;
        end
        t = toc; tproc = tproc + t; tic;
        fprintf('PostDark %gs. ',t)
    end
    % Final dark field
    if ~doDarkFieldCorrection
        darkVariant = 0;
    end
    switch darkVariant
        case 0
            darkMean = 0;
            fprintf('\n DARK FIELD IS SET TO ZERO.\n PROCESSING: ')
        case 1
            darkMean = preDarkMean;
            fprintf('\n ONLY PREDARK FIELD IS USED.\n PROCESSING: ')
        case 2
            darkMean = postDarkMean;
            fprintf('\n ONLY POSTDARK FIELD IS USED.\n PROCESSING: ')
        case 3
            darkMean = (preDarkMean + postDarkMean)/2;
    end
    %% PREFLAT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    flatVariant = 0;
    if NumPreFlat > 2
        NumPreFlat = mod(NumPreFlat,maxNumDarkFlat);
        stack = zeros(dim1,dim2,NumPreFlat-offset,'single');
        hpfPreFlat = zeros(1,size(stack,3),'single');
        parfor nn = 1:NumPreFlat-offset
            filename = sprintf('%s%s',InputPath,preFlatStruct(nn+offset).name);
            [stack(:,:,nn) hpfPreFlat(nn)] = FilterHotPixel(single(imread(filename,'tif','PixelRegion',PixelRegion)),hptFlat);
        end
        % Fourier space filtering
        if FSfiltering(2)
            stack = DataProc32ID_FilterFS(stack,'PreFlats');
        end
        preFlatMean = mean(stack,3);
        t = toc; tproc = tproc + t; tic;
        flatVariant = flatVariant + 1;
        fprintf('PreFlat %gs. ',t);
    end
    %% POSTFLAT
    if NumPostFlat > 2
        NumPostFlat = mod(NumPostFlat,maxNumDarkFlat);
        stack = zeros(dim1,dim2,NumPostFlat-offset,'single');
        hpfPostFlat = zeros(1,size(stack,3),'single');
        parfor nn = 1:NumPostFlat-offset
            filename = sprintf('%s%s',InputPath,postFlatStruct(nn+offset).name);
            [stack(:,:,nn) hpfPostFlat(nn)] = FilterHotPixel(single(imread(filename,'tif','PixelRegion',PixelRegion)),hptFlat);
        end
        % Fourier space filtering
        if FSfiltering(2)
            stack = DataProc32ID_FilterFS(stack,'PostFlats');
        end
        postFlatMean = mean(stack,3);
        t = toc; tproc = tproc + t; tic;
        flatVariant = flatVariant + 2;
        fprintf('PostFlat %gs. ',t);
    end
    switch flatVariant
        case 0
            fprintf('\n !!NO FLAT FIELDS FOUND. EXITING LOOP!!\n\n')
            return
        case 1
            flatMean = preFlatMean;
            fprintf('\n !ONLY PREFLATS FOUND! \n PROCESSING: ')
        case 2
            flatMean = postFlatMean;
            fprintf('\n !ONLY POSTFLATS FOUND! \n PROCESSING: ')
        case 3
            flatMean = (preFlatMean + postFlatMean)/2;
    end
    %% Modify dark field
    if doDarkFieldCorrection
        darkMeanAv = mean(darkMean(:));
        mask = ones(dim1,dim2);
        mask(flatMean < 2*mean(darkMeanAv)) = 0;
        mask = medfilt2(mask,[3 3],'symmetric');
        %mask = imfilter(mask,fspecial('gaussian',[3 3],10),'symmetric');
        mask = imfilter(mask,fspecial('disk',10),'symmetric');
        darkMean = darkMean.*mask;
        flatMean = flatMean - darkMean;
    end
    if doPreFilteringHorizontal
        flatMean = flatMean./repmat(mean(flatMean,2),[1 dim2]);
    end
    %% ROTATION AXIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Read first projection
    filename = sprintf('%s%s',InputPath,projStruct(1+offset).name);
    im = (FilterHotPixel(single(imread(filename,'tif','PixelRegion',PixelRegion)),hptProj) - darkMean)./flatMean;
    % Read last projection and flip
    filename = sprintf('%s%s',InputPath,projStruct(NumProj).name);
    im2 = (FilterHotPixel(single(imread(filename,'tif','PixelRegion',PixelRegion)),hptProj) - darkMean)./flatMean;
    im2 = fliplr(im2);
    % Crop vertically
    x = floor(3/10*dim1):ceil(dim1*7/10);
    im  = im(x,:);
    im2 = im2(x,:);
    % Compute axis
    out        = ImageCorrelation(im,im2,0,0);
    RotAxisPos = out.VerticalRotationAxisPosition;
    %% Make par file
    ParInputFilePrefix = sprintf('%sphase_',OutputPath);
    ParOutputFilePrefix = sprintf('%stomo%02u.vol',VolOutputPath,tomoInd);
    MakeParFile32ID(ParInputFilePrefix,ParOutputFilePrefix,[dim1 dim2],RotAxisPos,'NumberOfProjections',dim3,'OverallAngle',180*(1-offset/NumProj));
    fprintf('\n ROTATION AXIS POSITION: %g\n',RotAxisPos);
    %% PROJECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read, hot-pixel filter, flat field correction, phase retrieval
    stack = zeros(dim1,dim2,dim3,'single');
    hpfProj = zeros(1,dim3,'single');
    if useGPU
        darkMean = gpuArray(darkMean);
        flatMean = gpuArray(flatMean);
        phaseFilter = gpuArray(phaseFilter);
        parfor nn = 1:dim3
            filename = sprintf('%s%s',InputPath,projStruct(nn+offset).name);
            [im hpfProj(nn)] = FilterHotPixel(single(imread(filename,'tif','PixelRegion',PixelRegion)),hptProj);
            im = gpuArray(im);
            im = im - darkMean;
            if doPreFilteringHorizontal
                im = im./repmat(mean(im,2),[1 dim2]);
            end
            im = im./flatMean;
            if doPostFilteringHorizontal
                im = im./repmat(mean(im,2),[1 dim2]);
            end
            if doPostFilteringVertical
                im = im./repmat(mean(im,1),[dim1 1]);
            end
            if doMeanSubstraction
                im = im - mean(mean(im));
            end
            im = fft2(im);
            im = real(ifft2(im.*phaseFilter));
            edfwrite(sprintf('%sphase_%04u.edf',OutputPath,nn),gather(im)','float32');
        end
    else
        parfor nn = 1:dim3
            filename = sprintf('%s%s',InputPath,projStruct(nn+offset).name);
            [stack(:,:,nn) hpfProj(nn)] = FilterHotPixel(single(imread(filename,'tif','PixelRegion',PixelRegion)),hptProj);
            stack(:,:,nn) = stack(:,:,nn) - darkMean;
            if doPreFilteringHorizontal
                stack(:,:,nn) = stack(:,:,nn)./repmat(mean(stack(:,:,nn),2),[1 dim2]);
            end
            stack(:,:,nn) = stack(:,:,nn)./flatMean;
            if doPostFilteringHorizontal
                stack(:,:,nn) = stack(:,:,nn)./repmat(mean(stack(:,:,nn),2),[1 dim2]);
            end
            if doPostFilteringVertical
                 stack(:,:,nn) = stack(:,:,nn)./repmat(mean(stack(:,:,nn),1),[dim1 1]);
            end
            if doMeanSubstraction
                stack(:,:,nn) = stack(:,:,nn) - mean(mean(stack(:,:,nn)));
            end
        end
        stack = fft(stack,[],2);
        if doFilterSino
            stack = fft(stack,[],3);            
            stackz = median(stack(:,:,[1:3 end-1:end]),3);
            stack(:,:,1) = stackz;
            stack = ifft(stack,[],3);
        end
        stack = fft(stack,[],1);
        % Phase retrieval
        stack = stack.*repmat(phaseFilter,[1 1 dim3]);
        stack = ifft(ifft(stack,[],1),[],2);
        stack = real(stack);
        % Normalize slices        
        if doNormalizeSlices
            stack = stack./repmat(abs(mean(mean(stack,3),2)),[1 dim2 dim3]);
        end
        %% Save processed images as edf or tif
        switch writeFormat
            case 'edf'     
                parfor nn = 1:dim3
                    edfwrite(sprintf('%sphase_%04u.edf',OutputPath,nn),squeeze(stack(:,:,nn))','float32');
                end
            case {'tif','tiff'}
                parfor nn= 1:dim3
                    write32bitTIFfromSingle(sprintf('%sphase_%04u.edf',OutputPath,nn),im);
                end
        end
    end
    t = toc; tproc = tproc + t; tic;
    fprintf(' Sample & phase %gs. ',t); 
    %% Hot pixel statistics
    Numhpf = 0;
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
ttotal = ttotal + tpretomoloop + tpreprocloop + tproc;
fprintf('TOTAL PROCESSING TIME: %gs (PreTomoLoop %.2gs, PreProcLoop %.2gs, TomoProc: %gs)\n',ttotal,tpretomoloop,tpreprocloop,tproc);
%% Hot pixel statistics
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
    save([ParentOutputPath 'hpf.mat'],'hpf*Ar');
end
if PipeGPUtag && useGPU
    matlabpool close force local
end
clear all
end
%% END OF PRORGRAMME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
