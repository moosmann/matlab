function RecoLoopLCIstack(TomoSetsToProcess,doRingFilter,varargin)
% Main script for phase retrieval based on the script 'RecoLoop', but
% altered for 2-BM@APS beam time data. Calling script is intended to be
% RecoLoopAPS.
tic; telapsed = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Default arguments
if nargin < 1
    % Scalar (0,default), or 1x?-Vector Data sets to process. Default (0)
    % is processes all data sets found under 'IntParentPath', else
    % scalar or 1x?-Vector of the data set number(s) that should be processed.
    TomoSetsToProcess = 1;%input(sprintf('Tomo sets to process (row vector): '));
end
if nargin < 2
    doRingFilter = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input arguments
for nn = 1:2:numel(varargin)
    evalc([varargin{nn} ' = ' mat2str(varargin{nn+1})]);
    %assignin('caller',varargin{nn},varargin{nn+1});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('ParentPath','var')
    ParentPath = '/mnt/tomoraid-LSDF/tomo/APS_2BM_LifeCellImaging_GUP28266/lifecellimaging/';
end
if ~exist('DataSet','var')
    %DataSet = 'wildtype_30keV_05min_deadtime_08tomo_stage18p0_upwards_620mm_075ms';
    %DataSet = 'wildtype_30keV_10min_deadtime_25tomo_stage11p0_upwards_620mm_015ms';  
    DataSet = 'wildtype_30keV_10min_deadtime_20tomo_stage14p0_upwards_620mm_025ms';
end
if ~exist('PreProcessingFolder','var')
    %PreProcessingFolder = 'int_filtLineSectionMedWidthH063V001_noCropping_filtSino';
    PreProcessingFolder = 'int_filtLineSectionMedWidthH063V001_noCropping';
    %PreProcessingFolder = '';
end
if ~exist('EnergyDistancePixelsize','var')
    EnergyDistancePixelsize = [30 0.620 2.2e-6];
end
if ~exist('PhaseMethod','var')
    PhaseMethod = 'tie';
    %PhaseMethod = 'pctf';
    %PhaseMethod = 'pctfHalfSine';
end
if ~exist('RegPar','var')
    RegPar = 2.5;
end
if ~exist('BinaryFilterThreshold','var')
    BinaryFilterThreshold = 0.2;
end
if ~exist('Padding','var')
    % padding method is symmetric
    % 0: no padding, 1: pad to the next power of 2, 2: double padding
    % 2 avoid jumps at image boundaries and thus artifacts in FFT
    Padding = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dataType = 'single';
dataType = 'double';
% Experimental parameters.
Energy    = EnergyDistancePixelsize(1);
Distance  = EnergyDistancePixelsize(2);
Pixelsize = EnergyDistancePixelsize(3);
lambda    = EnergyConverter(Energy);
ArgPrefac = 2*pi*lambda*Distance/Pixelsize^2;
%% Folders and directories
CheckTrailingSlash(ParentPath)
IntParentPath = [ParentPath 'int/' DataSet '/' PreProcessingFolder];
CheckTrailingSlash(IntParentPath)
InputFolderNames = dir([IntParentPath 'tomo*']);
ParentPath = fileparts(ParentPath);
ParentOutputPath = [ParentPath '/phase/' DataSet '/' PreProcessingFolder];
CheckTrailingSlash(ParentOutputPath)
% Print info.
fprintf('\nPHASE RETRIEVAL\nDATA SET: ''%s''',DataSet)
fprintf('\nPATH TO TOMOGRAMS: %s',IntParentPath)
fprintf('\nTOMOGRAMS TO PROCESS:')
if TomoSetsToProcess == 0
    fprintf(' ALL\n')
else
    fprintf(' %s\n',mat2str(TomoSetsToProcess))
end
%% Loop over different data sets (frog embryo stages).
if TomoSetsToProcess == 0
    TomoSetsToProcess = 1:numel(InputFolderNames);
end
FolderPostfix = '';
if doRingFilter
    FolderPostfix = ['_filtRing' FolderPostfix];
end
%% Loop over projections
InputPath = [IntParentPath InputFolderNames(TomoSetsToProcess(1)).name '/'];
%% Read first image, set parameters and create the filter needed for phase retrieval.
im = pmedf_read_jm(sprintf('%sint_%04u.edf',InputPath,1),dataType);
[dim1,dim2]   = size(im); 
% Padding
switch lower(Padding)
    case {0,'no','nopadding','padno'}
        dimx          = dim1;
        dimy          = dim2;
        PaddingMethod = '_padNo';
    case {1,'npo2','nextpowerof2','padnpo2','padnextpowerof2'}
        dimx          = 2^nextpow2(dim1);
        dimy          = 2^nextpow2(dim2);
        PaddingMethod = '_padNPO2';
    case {2,'double','paddouble'}
        dimx          = 2*dim1;
        dimy          = 2*dim2;
        PaddingMethod = '_padDouble';
end
FolderPostfix = [PaddingMethod FolderPostfix];
xcut = 1+ceil((dimx-dim1)/2):floor((dimx+dim1)/2);
ycut = 1+ceil((dimy-dim2)/2):floor((dimy+dim2)/2);
% Preallocation
InputNameStruct = dir([InputPath 'int_*.edf']);
dim3 = numel(InputNameStruct);
stack = zeros(dimx,dimy,dim3,dataType);
[phaseFilter PostfixBinFilt] = PhaseFilter(PhaseMethod,[dimx dimy],EnergyDistancePixelsize,RegPar,BinaryFilterThreshold,dataType);
phaseFilter = repmat(phaseFilter,[1 1 dim3]);
fprintf(' INPUT PATH:  %s\n',InputPath)
fprintf(' PARENT OUTPUT PATH: %s\n',ParentOutputPath)
fprintf(' SIZE OF PROJECTION: hor x vert = %u x %u (unpadded), %u x %u (padded)\n',dim1,dim2,dimx,dimy)
fprintf(' SIZE OF STACK: dim1 x dim2 x dim3 = %u x %u x %u\n',size(stack))
fprintf(' REGULARIZATI0N PARAMETER: alpha = %g, 1/2/( 2*pi*lambda*z/dx^2 * xi^2/2 + 10^-alpha ) = 1/2/( %g * [-1/2,-1/2-1/N] + %g )\n',RegPar,ArgPrefac,10^-RegPar)
%% Phase retrieval filters
for kk = TomoSetsToProcess
    InputPath = [IntParentPath InputFolderNames(kk).name '/'];
    %% Create folders for retrieved-phase images will be stored.
    OutputPath = sprintf('%s%s_alpha%3.2f%s%s/%s/',ParentOutputPath,PhaseMethod,RegPar,PostfixBinFilt,FolderPostfix,InputFolderNames(kk).name);
    OutputPath = regexprep(OutputPath,'\.','p');
    CheckAndMakePath(OutputPath)
    %% Start looping.
    fprintf('\n PROCESSING TOMO SET: %s\n',InputFolderNames(kk).name)
    fprintf(' OUTPUT PATH: %s\n',OutputPath)
    %% Read projections into stack
    telapsed = toc; tic
    parfor nn = 1:dim3
        %% Reading, padding and FT of data.
        %stack(:,:,nn) = pmedfread([InputPath InputNameStruct(nn).name],dataType)
        stack(:,:,nn) = padarray(padarray((pmedf_read_jm([InputPath InputNameStruct(nn).name],dataType)),[ceil((dimx-dim1)/2),ceil((dimy-dim2)/2)],'symmetric','pre'),[floor((dimx-dim1)/2),floor((dimy-dim2)/2)],'symmetric','post');
    end
    tread = toc; telapsed = telapsed + tread; tic
    %% Filter stack
    stack = fft(stack,[],1);
    if doRingFilter
        stack = fft(stack,[],3);
        stackz = median(stack(:,:,[1:3 end-1:end]),3);
        stack(:,:,1) = stackz;
        stack = ifft(stack,[],3);
    end
    stack = fft(stack,[],2);
    stack = stack.*phaseFilter;
    %stack = real(ifft(ifft(stack,[],2,'symmetric'),[],1,'symmetric'));
    stack = real(ifft(ifft(stack,[],2),[],1));
    tfilter = toc; telapsed = telapsed + tfilter; tic
    parfor nn = 1:dim3
        %% SAVE
        edfwrite(sprintf('%sphase_%04u.edf',OutputPath,nn),squeeze(stack(xcut,ycut,nn)),'float32');
    end
    tsave = toc; telapsed = telapsed + tsave; tic
   fprintf('Processing time: Read %g s. Filter %g s. Save %g s.\n',tread,tfilter,tsave) 
end
telapsed = telapsed + toc;
fprintf('Total processing time: %g min, %g h\n',telapsed/60,telapsed/3600)
