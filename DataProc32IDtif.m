function DataProc32IDtif(TomoSetsToProcess,useGPU,showFigures,varargin)
% Data preprocessing for the experiment 'Comparative Life Cell Imaging' at
% 32-ID-C at APS in October 2012

tic; tpretomoloop = 0; tpreprocloop = 0; tproc = 0; ttotal = 0;
%% Default arguments
if nargin < 1
    % Index of tomogram as it appears in the filename. It's not the
    % absolute number of tomograms
    TomoSetsToProcess = 1;%{};%:13
end
if nargin < 2
    useGPU = 1;
end
if nargin < 3
    showFigures = 0;
end
%% Input arguments
for nn = 1:2:numel(varargin)
    evalc([varargin{nn} ' = ' mat2str(varargin{nn+1})]);
    %assignin('caller',varargin{nn},varargin{nn+1});
end
%% Set default parameter
if ~exist('ParentPath','var')
    ParentPath = '/mnt/tomoraid-LSDF/tomo/APS_32ID-C_LifeCellImaging_GUP31523_2012-10-13/savedLocally';
    %ParentPath = '/mnt/tomoraid-LSDF/tomo/APS_32ID-C_LifeCellImaging_GUP31523_2012-10-13/savedExternally/';
    %ParentPath = '/mnt/gpuNodeStorage/scratch/moosmann/APS_32ID-C_LifeCellImaging_GUP31523_2012-10-13/savedLocally');
end
if ~exist('DataSet','var')
    %DataSet = 'Oct11_00-30_wildtype_stage11p0_34p5keV_0700mm_15ms_0834proj_scantime50s_deadtime08min_20ms_open_40ms_close';
    %DataSet = 'Oct10_21-48_wildtype_stage11p0_34p5keV_0700mm_20ms_0834proj_scantime50s_deadtime10min_24ms_open_36ms_close';
    %DataSet = 'Oct12_04-52_wildtype_stage12p0_34p5keV_0700mm_30ms_0834proj_scantime50s_deadtime10min_33ms_open_27ms_close';
    DataSet = 'Oct11_02-48_wildtype_stage11p5_34p5keV_0700mm_20ms_0834proj_scantime50s_deadtime10p5min_20ms_open_40ms_close';
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
    PixelRegion = {};%{[400 1200],[61 1400]};
    %PixelRegion = {[201 1274],[1 1488]};
end
if ~exist('EnergyDistancePixelsize','var')
    EnergyDistancePixelsize = [34.5 0.7 1.1e-6];
end
if ~exist('RegPar','var')
    RegPar = 2.5;
end
% if ~exist('TomoSetsToProcess','var')
%     TomoSetsToProcess = 3;%:13;
% end
%% Open pool of MATLAB sessions for parallel computation
if matlabpool('size') == 0
    matlabpool('open')
end
%% Parameter
hptDark = HotPixelThreshold(1);
hptFlat = HotPixelThreshold(2);
hptProj = HotPixelThreshold(3);
offset = 1; % Dismiss first image of each sequence due to camera start up time
Energy    = EnergyDistancePixelsize(1);
Distance  = EnergyDistancePixelsize(2);
Pixelsize = EnergyDistancePixelsize(3);
lambda    = 6.62606896e-34*299792458/(Energy*1.60217733e-16);
% Prefactor needed for TIE and CTF retrieval.
ArgPrefac = 2*pi*lambda*Distance/Pixelsize^2;
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
fprintf('\nSTARTING DATA PROCESSING\nDATA SET: %s\n',DataSet(1:end-1))
fprintf(' DATA PATH: %s\n',DataPath)
%fprintf('INPUT PATH: %s\n',InputPath)
% Number of tomograms found
TomoStruct = dir(sprintf('%s*postDark_00019.tif',InputPath));
NumTomosFound = numel(TomoStruct);
fprintf(' TOMOGRAMS FOUND: %u\n',NumTomosFound)
fprintf(' HOT PIXEL FILTER THRESHOLOD: : [dark flat proj] = [%g %g %g]\n',HotPixelThreshold)
%% Loop over tomograms
if isempty(TomoSetsToProcess)
    for nn = 1:NumTomosFound
        TomoSetsToProcess(nn) = str2double(TomoStruct(nn).name(6));
    end
    TomoSetsToProcess = sort(TomoSetsToProcess);
end
NumTomos = numel(TomoSetsToProcess);
% Read image dimensions for preallocation
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
fprintf(' IMAGE DIMENSIONS: [ver x hor] = [%u x %u]\n',dim1,dim2)
fprintf(' REGULARIZATI0N PARAMETER: 10^-%g = %g',RegPar,10^-RegPar)
%% Phase retrieval filters
% Fourier coordinates.
xi  = single(-1/2:1/dim2:1/2-1/dim2);
eta = single(-1/2:1/dim1:1/2-1/dim1);
[xiquad sinxiquad]   = meshgrid(xi,eta);
xiquad = fftshift(xiquad.^2 +sinxiquad.^2);
% sinxiquad  = sin(ArgPrefac*xiquad/2);
if useGPU
    tie = gpuArray(1./(ArgPrefac*xiquad + 10^-RegPar));
else
    tie = 1./(ArgPrefac*xiquad + 10^-RegPar);
    % InverseSine   = 1./(2*sign(sinxiquad).*(abs(sinxiquad))+10^-RegPar);
end
% Preallocation
hpfPreDarkAr = cell(NumTomos,1);
hpfPreFlatAr = cell(NumTomos,1);
hpfProjAr    = cell(NumTomos,1);
tpretomoloop = tpretomoloop + toc;tic
%% Start looping
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
    if NumPreDark*NumPreFlat*NumProj*NumPostFlat*NumPostFlat == 0
        fprintf('\n\n !LOOP OVER TOMOGRAMS STOPPED! (TOMOGRAM NUMBER: %u. LOOP INDEX: %u.)\n\n',tomoInd,tt)
        return
    end
    %% Paths
    PhaseFolderName = sprintf('FDcor_tie_regPar%3.2f/',RegPar);
    PhaseFolderName = regexprep(PhaseFolderName,'\.','p');
    OutputFolder = sprintf('%stomo%02u/',PhaseFolderName,tomoInd);
    fprintf(' OUTPUT FOLDER: %s\n',OutputFolder);
    OutputPath = sprintf('%s%s',ParentOutputPath,OutputFolder);
    CheckAndMakePath(OutputPath);
    VolOutputPath = sprintf('%s%s',ParentVolPath,PhaseFolderName);
    CheckAndMakePath(VolOutputPath);
    %% preDarks: read, hot-pixel filter
    tpreprocloop = tpreprocloop + toc; tic
    stack = zeros(dim1,dim2,NumPreDark-offset,'single');
    hpfPreDark = zeros(1,size(stack,3),'single');
    parfor nn = 1:NumPreDark-offset
        filename = sprintf('%s%s',InputPath,preDarkStruct(nn+offset).name);
        [stack(:,:,nn) hpfPreDark(nn)] = FilterHotPixel(single(imread(filename,'tif','PixelRegion',PixelRegion)),hptDark);
    end
    darkMean = mean(stack,3);
    t = toc; tproc = tproc + t; tic;
    fprintf('Made PreDark: %g s. ',t)
    %% preFlats: read, hot-pixel filter
    stack = zeros(dim1,dim2,NumPreFlat-offset,'single');
    %     mask = ones(dim2,size(stack,3),'single');
    %     mask(2:end,1) = 0; mask(1,2:end) = 0;
    hpfPreFlat = zeros(1,size(stack,3),'single');
    parfor nn = 1:NumPreFlat-offset
        filename = sprintf('%s%s',InputPath,preFlatStruct(nn+offset).name);
        [stack(:,:,nn) hpfPreFlat(nn)] = FilterHotPixel(single(imread(filename,'tif','PixelRegion',PixelRegion)),hptFlat);
    end
    %% preFlats: Fourier space filtering
    if FSfiltering(2)
        parfor nn = 1:dim1
            im = squeeze(stack(nn,:,:));
            im = fft2(im);
            im(2:end,1) = 0; im(1,2:end) = 0;
            stack(nn,:,:) = real(ifft2(im));
            %stack(nn,:,:) = real(ifft2(mask.*fft2(squeeze(stack(nn,:,:)))));
        end
        fprintf('FS Filtered PreFlats. ')
    end
    flatMean = mean(stack,3) -darkMean;
    t = toc; tproc = tproc + t; tic;
    fprintf('Made PreFlat: %g s. ',t);
    %% postDarks: read, hot-pixel filter
    %     parfor nn = 1+offset:NumPostDark
    %         filename = sprintf('%s%s',InputPath,postDarkStruct(nn).name);
    %         [stack(:,:,nn-offset) hpfPostDark] = FilterHotPixel(single(imread(filename,'tif')),hptDark);
    %     end
    %     darkMean = mean(stack,3);
    %   stack = zeros(dim1,dim2,NumPostFlat-offset,'single');
    %% postFlats: read, hot-pixel filter
    %     parfor nn = 1+offset:NumPostFlat
    %         filename = sprintf('%s%s',InputPath,postFlatStruct(nn).name);
    %         [stack(:,:,nn-offset) hpfPostFlat] = FilterHotPixel(single(imread(filename,'tif')),hptFlat);
    %     end
    %% proj: read, hot-pixel filter
    %NumProj = floor(NumProj/NumProj*NumPreFlat);
    %     stack = zeros(dim1,dim2,NumProj-offset,'single');
    %     hpfProj = zeros(1,NumProj-offset,'single');
    %     parfor nn = 1:NumProj-offset
    %         filename = sprintf('%s%s',InputPath,projStruct(nn+offset).name);
    %         [stack(:,:,nn) hpfProj(nn)] = FilterHotPixel(single(imread(filename,'tif','PixelRegion',PixelRegion)),hptProj);
    %         stack(:,:,nn) = (stack(:,:,nn) - darkMean)./flatMean;
    %     end
    %     parfor nn = 1:NumProj-offset
    %         stack(:,:,nn) = stack(:,:,nn)./flatMean;
    %     end
    %% Loop Projections: Fourier space filtering
    if FSfiltering(3)
        parfor nn = 1:dim1
            im = squeeze(stack(nn,:,:));
            im = fft2(im);
            im(1:end,1) = 0; im(1,1:end) = 0;
            stack(nn,:,:) = real(ifft2(im));
            %stack(nn,:,:) = real(ifft2(mask.*fft2(squeeze(stack(nn,:,:)))));
        end
        fprintf('FS Filtered Projs. ')
    end
    %     t = toc; tproc = tproc + t; tic
    %     fprintf('Flat-and-Darkfield correction: %g s. ',t)
    %% proj: read, hot-pixel filter, flat field correction, phase retrieval
    % sequential processing
    %clear stack im im2;
    hpfProj = zeros(1,NumProj-offset,'single');
    if useGPU
        parfor nn = 1:NumProj-offset
            %fprintf('%4u',nn)
            filename = sprintf('%s%s',InputPath,projStruct(nn+offset).name);
            %im = FilterHotPixel(single(imread(filename,'tif','PixelRegion',PixelRegion)),hptProj);
            [im hpfProj(nn)] = FilterHotPixel(single(imread(filename,'tif','PixelRegion',PixelRegion)),hptProj);
            im = gpuArray(im);
            im = (im - darkMean)./flatMean;
            im = fft2(im - mean(im(:)));
            im = gather(real(ifft2(im.*tie)));
            write32bitTIFfromSingle(sprintf('%sphase_%04u.tif',OutputPath,nn),im);
            %edfwrite(sprintf('%sphase_%04u.edf',OutputPath,nn),gather(im)','float32');
            
        end
    else
        parfor nn = 1:NumProj-offset
            %fprintf('%4u',nn)
            filename = sprintf('%s%s',InputPath,projStruct(nn+offset).name);
            %im = FilterHotPixel(single(imread(filename,'tif','PixelRegion',PixelRegion)),hptProj);
            [im hpfProj(nn)] = FilterHotPixel(single(imread(filename,'tif','PixelRegion',PixelRegion)),hptProj);
            im = (im - darkMean)./flatMean;
            im = fft2(im - mean(im(:)));
            im = real(ifft2(im.*tie));
            write32bitTIFfromSingle(sprintf('%sphase_%04u.tif',OutputPath,nn),im);
            %edfwrite(sprintf('%sphase_%04u.edf',OutputPath,nn),im','float32');
        end
    end
    % 3D processing
    %     if useGPU
    %         parfor nn = 1:NumProj-offset
    %             fprintf('%4u',nn)
    %             im = gpuArray(squeeze(stack(:,:,nn)));
    %             im = fft2(im - mean(im(:)));
    %             im = real(ifft2(im.*tie));
    %             %stack(:,:,nn) = im;
    %             edfwrite(sprintf('%sphase_%04u.edf',OutputPath,nn),gather(im)','float32');
    %         end
    %     else
    %         parfor nn = 1:NumProj-offset
    %             fprintf('%4u',nn)
    %             im = squeeze(stack(:,:,nn));
    %             im = fft2(im - mean(im(:)));
    %             im = real(ifft2(im.*tie));
    %             %stack(:,:,nn) = im;
    %             edfwrite(sprintf('%sphase_%04u.edf',OutputPath,nn),im','float32');
    %         end
    %     end
    t = toc; tproc = tproc + t; tic;
    fprintf(' Projections (preprocessing & phase retrieval: %g s. ',t);
    %% Rotation axis
    % Read first projection
    filename = sprintf('%s%s',InputPath,projStruct(1+offset).name);
    im = (FilterHotPixel(single(imread(filename,'tif','PixelRegion',PixelRegion)),hptProj) - darkMean)./flatMean;
    % Read last projection
    filename = sprintf('%s%s',InputPath,projStruct(NumProj).name);
    im2 = (FilterHotPixel(single(imread(filename,'tif','PixelRegion',PixelRegion)),hptProj) - darkMean)./flatMean;
    im2 = fliplr(im2);
    % Crop vertically
    x = floor(2/10*dim1):ceil(dim1*9/10);
    im  = im(x,:);
    im2 = im2(x,:);
    % Compute axis
    out         = ImageCorrelation(im,im2,0,0);
    RotAxisPos  = out.VerticalRotationAxisPosition;
    %% Make par file
    ParInputFilePrefix = sprintf('%sphase_',OutputPath);
    ParOutputFilePrefix = sprintf('%stomo%02u.vol',VolOutputPath,tomoInd);
    MakeParFile32ID(ParInputFilePrefix,ParOutputFilePrefix,[dim1 dim2],RotAxisPos,'NumberOfProjections',NumProj-offset);
    t = toc; tproc = tproc + t; tic;
    fprintf(' Rotation axis & par file: %g s. ',t);
    hpfPreDarkAr{tt} = hpfPreDark;
    hpfPreFlatAr{tt} = hpfPreFlat;
    hpfProjAr{tt}    = hpfProj;
    fprintf('\n ROTATION AXIS POSITION; %g\n',RotAxisPos);
end
ttotal = ttotal + tpretomoloop + tpreprocloop + tproc;
fprintf('TOTAL PROCESSING TIME in s: %g (PreTomoLoop %.2g, PreProcLoop %.2g, TomoProc: %g)\n',ttotal,tpretomoloop,tpreprocloop,tproc);
%% For analysis issues
if showFigures
    figure('Name','Percentage of filtered hot pixels')
    subplot(3,1,1), plot(100*[hpfPreDarkAr{:}],'x'), title('PreDark')
    subplot(3,1,2), plot(100*[hpfPreFlatAr{:}],'x'), title('PreFlat')
    subplot(3,1,3), plot(100*[hpfProjAr{:}],'x'), title('Projections')
end
if exist('hpf','var')
    save([ParentOutputPath 'hpf.mat'],'hpfPreDarkAr','hpfPreFlatAr','hpfProjAr');
end
end
