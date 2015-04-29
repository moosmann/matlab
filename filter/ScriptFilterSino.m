function [sino,sino_unfilt] = ScriptFilterSino(RowsToRead,InputDir,ImageNamePattern,OutputDirSino,OutputDirImages,ProjectionsToRead,HotPixThresh)
% Make sinogram

%% Default arguments
if nargin < 1
    %vector containing the row number that should be read
    RowsToRead = [1 40 94 138 176 223 270 324 404 453 519 586 664 720 773 860 943 994]; 
end
if nargin < 2
    %string of the absolute or relative input path. Default is the current folder
    %InputDir = '/mnt/tomoraid-LSDF/tomo/APS_2BM_LifeCellImaging_GUP28266/lifecellimaging/test/data/wildtype_30keV_10min_deadtime_25tomo_stage11p0_upwards_620mm_015ms';
    InputDir = '/mnt/tomoraid-LSDF/tomo/APS_2BM_LifeCellImaging_GUP28266/lifecellimaging/int/wildtype_30keV_10min_deadtime_25tomo_stage11p0_upwards_620mm_015ms/int_filtLineSectionMedWidthH063V001_noCropping/tomo01';
end
if nargin < 3
    %string containg the naming pattern
    ImageNamePattern = []; 
end
if nargin < 4
    % String of the output path where the sinograms will be stored
    OutputDirSino = []; 
end
if nargin < 5
    OutputDirImages = [];
end
if nargin < 6
    ProjectionsToRead = 1:1200;
end
if nargin < 7
    HotPixThresh = 0;
end
% if nargin < 7
%     CircularFilterWidth = 2;
% end
% if nargin < 8
%     CrossFilterWidth = 2;
% end
% if nargin < 9
%     LineFilterWidth = 2;
% end
% if nargin < 10
%     theta = 76.15;
% end
% if nargin < 11
%     DoFFTshift = 1;
% end
%% Body
tic;
% Check InputDir string
if InputDir(end) ~= '/'
   InputDir = [InputDir '/'];
end
% Read image struct
if isempty(ImageNamePattern)
    ImStruct = dir([InputDir '*.*']);
    ImStruct(1:2) = [];
else
    ImStruct = dir([InputDir ImageNamePattern]);
end
NumIm = numel(ImStruct);
% Check image file format
ImFormat = ImStruct(1).name(end-2:end);
% Get image dimensions
if strcmp(ImFormat,'edf')
    im = edfread([InputDir ImStruct(1).name])';
else
    im = double(imread([InputDir ImStruct(1).name]));
end
[dimx,dimy] = size(im);
%clear im;
% Assing default row for sinogram
if RowsToRead(1) == 0
    RowFirst = ceil(dimx/2);
    RowLast  = RowFirst;
else
    RowFirst = RowsToRead(1);
    RowLast  = RowsToRead(end);
end
% Decleare images that will be used for the tomogram
if isempty(ProjectionsToRead)
    ProjectionsToRead = 1:1:NumIm;
end
NumProjs = numel(ProjectionsToRead);
FileNames = {ImStruct(ProjectionsToRead).name};
NumRows = numel(RowsToRead);
% Preallocate sinogram
sino = zeros(NumProjs,dimy,NumRows);
%% Loop over images and make sinogram
if strcmp(ImFormat,'edf')
    parfor nn = 1:NumProjs
        if HotPixThresh == HotPixThresh
            sino(nn,:,:) = edfread([InputDir FileNames{nn}],1:dimy,RowsToRead);
        else
            sino(nn,:,:) = FilterHotPixel(edfread([InputDir FileNames{nn}],1:dimy,RowsToRead),HotPixThresh);
        end
    end
else
    parfor nn = 1:NumProjs
        if HotPixThresh == HotPixThresh
            sino(nn,:,:) = double(imread([InputDir FileNames{nn}],'PixelRegion',{[RowFirst RowLast], [1 dimy]}))';
        else
            sino(nn,:,:) = FilterHotPixel(double(imread([InputDir FileNames{nn}],'PixelRegion',{[RowFirst RowLast], [1 dimy]}))',HotPixThresh);
        end
    end
end
% Print info
fprintf('Read rows %u to %u of %u %s-images of dimension [%u x %u] into sinogram of dimension [%u x %u x %u] in %fs.\n',...
    RowFirst,RowLast,NumProjs,ImFormat,dimx,dimy,size(sino,1),size(sino,2),size(sino,3),toc)
% %% Filter sinogram
% % Create mesh
% [x y] = meshgrid(-dimy/2:dimy/2-1,-dimx/2:dimx/2-1);
% %% Circular center-symmetric high-pass filter in Fourier space
% CircularFilter = 1-exp(-(x.^2+y.^2)/CircularFilterWidth^2/2);
% CrossFilter = (1-exp(-x.^2/CrossFilterWidth^2/2)).*(1-exp(-y.^2/CrossFilterWidth^2/2));
% LineFilter = ones(dimx,dimy);
% LineFilter(floor(abs(x-cotd(theta)*y))<LineFilterWidth) = 0;
% %LineFilter = fftshift(LineFilter);
% Filter = CrossFilter.*CircularFilter.*LineFilter;
% if DoFFTshift == 1
%     Filter     = fftshift(Filter);
% end

% % Filter image
% sinof = real(ifft2(Filter.*fft2(sino)));
tic
%% Normalize sino
%[dim1 , ~, dim3] = size(sino);
if nargout > 1
    sino_unfilt = sino;
end
parfor nn = 1:NumRows
    sino(:,:,nn) = sino(:,:,nn)./repmat(median(squeeze(sino(:,:,nn))),[NumProjs 1]);
end
fprintf('Time to normalize sinograms: %fs.\n',toc)
%% Save sinogram
if ~isempty(OutputDirSino)
    tic
    if OutputDirSino(end)~='/'
        OutputDirSino = [OutputDirSino '/'];
    end
    if ~exist(OutputDirSino,'dir')
        mkdir(OutputDirSino)
    end
    parfor nn = 1:NumRows
        edfwrite(sprintf('%ssino_%04u.edf',OutputDirSino,RowsToRead(nn)),sino(:,:,nn)','float32');
    end
    fprintf('Time to save sinograms: %fs.\n',toc)
end
%% Save images
if ~isempty(OutputDirImages)
    tic
    if OutputDirImages(end)~='/'
        OutputDirImages = [OutputDirImages '/'];
    end
    if ~exist(OutputDirImages,'dir')
        mkdir(OutputDirImages)
    end
    parfor nn = 1:NumProjs
        edfwrite(sprintf('%sint_%04u.edf',OutputDirImages,nn),squeeze(sino(nn,:,:)),'float32');
    end
    fprintf('Time to save projections: %fs.\n',toc)
end
