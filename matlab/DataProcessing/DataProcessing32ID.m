function DataProcessing32ID(varargin)
% Data preprocessing for the experiment 'Comparative Life Cell Imaging' at
% 32-ID-C at APS in October 2012

tic; tpre = 0; ttotal = 0; tpreloop = 0; tprocess = 0;
%% Open pool of MATLAB sessions for parallel computation
if matlabpool('size') == 0
    matlabpool('open')
end
%% Input arguments
for nn = 1:2:numel(varargin)
    evalc([varargin{nn} ' = ' mat2str(varargin{nn+1})]);
    %assignin('caller',varargin{nn},varargin{nn+1});
end
%% Set default parameter
if ~exist('ParentPath','var')
    ParentPath = '/mnt/tomoraid-LSDF/tomo/APS_32ID-C_LifeCellImaging_GUP31523_2012-10-13/savedLocally';
    %ParentPath = 'APS_32ID-C_LifeCellImaging_GUP31523_2012-10-13/savedLocally');
    %ParentPath = '/mnt/gpuNodeStorage/scratch/moosmann/APS_32ID-C_LifeCellImaging_GUP31523_2012-10-13/savedLocally');
end
if ~exist('DataSet','var')
    %DataSet = 'Oct11_00-30_wildtype_stage11p0_34p5keV_0700mm_15ms_0834proj_scantime50s_deadtime08min_20ms_open_40ms_close';
    DataSet = 'Oct12_11-45_wildtype_stage22p0_34p5keV_0700mm_30ms_0834proj_scantime50s_deadtime10min_33ms_open_27ms_close';
end
if ~exist('TomoSetsToProcess','var')
    TomoSetsToProcess = 2;
end
if ~exist('HotPixelThreshold','var')
    HotPixelThreshold = [1.2 1.1 1.1];
end
if ~exist('PixelRegion','var')
    %PixelRegion = {[400 1200],[61 1400]};
    PixelRegion = {};
end
% PixelRegion = {ROWS, COLS};
% ROWS = [START STOP]; COLS = [START STOP];
if ~exist('EnergyDistancePixelsize','var')
    EnergyDistancePixelsize = [34.5 0.7 1.1e-6];
end
if ~exist('RegPar','var')
    RegPar = 2.5;
end
if ~exist('GPU','var')
    GPU = true;
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
fprintf('\nSTARTING DATA PREPROCESSING\nDATA SET: %s\n',DataSet(1:end-1))
fprintf('DATA PATH: %s\n',DataPath)
%fprintf('INPUT PATH: %s\n',InputPath)
% Number of tomograms found
TomoStruct = dir(sprintf('%s*postDark_00019.tif',InputPath));
NumTomos = numel(TomoStruct);
fprintf('TOMOGRAMS FOUND: %u\n',NumTomos)
fprintf('HOT PIXEL FILTER THRESHOLOD: : [dark flat proj] = [%g %g %g]\n',HotPixelThreshold)
%% Loop over tomograms
if isempty(TomoSetsToProcess)
    TomoSetsToProcess = 1:NumTomos;
end
% Read image dimensions for preallocation
if isempty(PixelRegion)
    filename = sprintf('%sproj_0preDark_00000.tif',InputPath);
    iminfo = imfinfo(filename,'tif');
    dim1 = iminfo.Height;
    dim2 = iminfo.Width;
    PixelRegion = {[1 dim1],[1 dim2]};
else
    dim1 = PixelRegion{1}(2)-PixelRegion{1}(1)+1;
    dim2 = PixelRegion{2}(2)-PixelRegion{2}(1)+1;
end
%% Phase retrieval filters
% Fourier coordinates.
xi  = single(-1/2:1/dim2:1/2-1/dim2);
eta = (-1/2:1/dim1:1/2-1/dim1);
[xi eta]   = meshgrid(xi,eta);
xiquad = fftshift(xi.^2 +eta.^2);
tie = 1./(ArgPrefac*xiquad + 10^-RegPar);
%sinxiquad  = sin(ArgPrefac*xiquad/2);
%InverseSine   = 1./(2*sign(sinxiquad).*(abs(sinxiquad))+10^-RegPar);
fprintf('IMAGE DIMENSIONS: [ver x hor] = [%u x %u]\n',dim1,dim2)
fprintf('\n')
fprintf('REGULARIZATI0N PARAMETER: 10^-%g = %g\n',RegPar,10^-RegPar)
%%
tpre = tpre + toc;
% Start looping
for tomoInd = TomoSetsToProcess
    tic
    fprintf('\nPROCESSING TOMO INDEX: %u\n',tomoInd)
    % Check number of flats, darks, and projections
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
    fprintf('IMAGES FOUND: %u preDark, %u preFlat, %u proj, %u postFlat, %u postDark\n', ...
        NumPreDark,NumPreFlat,NumProj,NumPostFlat,NumPostDark)
    % Output path
    OutputPath = [ParentPath 'phase/'];
    PhaseMethod = sprintf('tie_alpha%3.2f',RegPar);
    PhaseMethod = regexprep(PhaseMethod,'\.','p');
    OutputPath = sprintf('%s%s%s/tomo%02u/',OutputPath,DataSet,PhaseMethod,tomoInd);
    if ~exist(OutputPath,'dir')
        mkdir(OutputPath);
    end
    fprintf('OUTPUT PATH: %s\n',OutputPath)
    cd(OutputPath)
    tpreloop = tpreloop + toc; tic
    %% Loop preDarks: read and hot-pixel filter
    fprintf('Processing: PreDarks. ')
    stack = zeros(dim1,dim2,NumPreDark-offset,'single');
    hpfPreDark = zeros(1,size(stack,3),'single');
    parfor nn = 1:NumPreDark-offset
        filename = sprintf('%s%s',InputPath,preDarkStruct(nn+offset).name);
        [stack(:,:,nn) hpfPreDark(nn)] = FilterHotPixel(single(imread(filename,'tif','PixelRegion',PixelRegion)),hptDark);
    end
    darkMean = mean(stack,3);
    %% Loop preFlats: read and hot-pixel filter
    fprintf('PreFlats. ')
    stack = zeros(dim1,dim2,NumPreFlat-offset,'single');
    hpfPreFlat = zeros(1,size(stack,3),'single');
    parfor nn = 1:NumPreFlat-offset
        filename = sprintf('%s%s',InputPath,preFlatStruct(nn+offset).name);
        [stack(:,:,nn) hpfPreFlat(nn)] = FilterHotPixel(single(imread(filename,'tif','PixelRegion',PixelRegion)),hptFlat);
        stack(:,:,nn) = stack(:,:,nn) - darkMean;
    end
    tprocess = tprocess + toc; tic
    flatMean = mean(stack,3);
    %% Loop proj: read and hot-pixel filter
    fprintf('Projections.')
    %NumProj = floor(NumProj/NumProj*NumPreFlat);
    hpfProj = zeros(1,NumProj-offset,'single');
    if GPU
        flatMean = gpuArray(flatMean);
        darkMean = gpuArray(darkMean);
        tie = gpuArray(tie);
        parfor nn = 1:NumProj-offset
            filename = sprintf('%s%s',InputPath,projStruct(nn+offset).name);
            [im hpfProj(nn)] = FilterHotPixel(single(imread(filename,'tif','PixelRegion',PixelRegion)),hptProj);
            im = gpuArray(im);
            im = (im - darkMean)./flatMean;
            filename = sprintf('%sint_%04u.tif',OutputPathInt,nn);
            write32bitTIF(filename,im')            
            im = fft2(im - mean(im(:)));
            im = real(ifft2(tie.*im));
            filename = sprintf('%sphase_%04u.tif',OutputPath,nn);
            write32bitTIF(filename,gather(im)')            
            %edfwrite(filename,gather(im'),'float32');
        end
    else
        parfor nn = 1:NumProj-offset
            filename = sprintf('%s%s',InputPath,projStruct(nn+offset).name);
            [im hpfProj(nn)] = FilterHotPixel(single(imread(filename,'tif','PixelRegion',PixelRegion)),hptProj);
            im = (im - darkMean)./flatMean;
            filename = sprintf('%sint_%04u.tif',OutputPathInt,nn);
            write32bitTIF(filename,im')
            im = fft2(im - mean(im(:)));
            im = real(ifft2(tie.*im));
            filename = sprintf('%sphase_%04u.tif',OutputPath,nn);
            write32bitTIF(filename,im')            
            %edfwrite(sprintf('%sphase_%04u.edf',OutputPath,nn),im','float32');
        end
    end
    tprocess = tprocess + toc; tic
end
ttotal = ttotal + tpre + tpreloop + tprocess;
fprintf('\nPreTomoLoop %.2g, PreProjLoop %.2g, ProjLoop: %g.\n',tpre,tpreloop,tprocess)
fprintf('TOTAL PROCESSING TIME: %g s\n',ttotal)
%% For analysis issues
figure('Name','Percentage of filtered hot pixels')
subplot(3,1,1), plot(100*hpfPreDark,'x'), title('PreDark')
subplot(3,1,2), plot(100*hpfPreFlat,'x'), title('PreFlat')
subplot(3,1,3), plot(100*hpfProj,'x'), title('Projections')
end

