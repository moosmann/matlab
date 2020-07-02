function [sino,ImStruct] = MakeSino(RowsToRead,InputDir,ImageNamePattern,OutputDir,OutputFormat,ProjectionsToRead,HotPixThresh)
% Make sinogram

%% Default arguments
if nargin < 1
    RowsToRead = 0; %vector containing the row number that should be read
end
if nargin < 2
    InputDir = './';%string of the absolute or relative input path. Default is the current folder
end
if nargin < 3
    ImageNamePattern = []; %string containg the naming pattern
end
if nargin < 4
    OutputDir = []; % String of the output path where the sinograms will be stored
end
if nargin < 5
    OutputFormat = 'tif';
end
if nargin < 6
    ProjectionsToRead = [];
end
if nargin < 7
    HotPixThresh = 0.00;
end
%% Body
tic;
% Check InputDir string
CheckTrailingSlash(InputDir);
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
[dimx, dimy] = size(im);
clear im;
% Assing default row for sinogram
if RowsToRead(1) == 0
    RowFirst = ceil(dimx/2);
    RowLast  = RowFirst;
else
    RowFirst = RowsToRead(1);
    RowLast  = RowsToRead(end);
end
% Declare images that will be used for the tomogram
if isempty(ProjectionsToRead)
    ProjectionsToRead = 1:1:NumIm;
end
NumImRead = numel(ProjectionsToRead);
FileNames = {ImStruct(ProjectionsToRead).name};
%% Loop over images and make sinogram
if strcmp(ImFormat,'edf')
    parfor nn = 1:NumImRead%:-1:1
        sino(nn,:,:) = FilterPixel(edfread([InputDir FileNames{nn}],1:dimy,RowsToRead),HotPixThresh);
    end
else
    parfor nn = 1:NumImRead%:-1:1
        sino(nn,:,:) = FilterPixel(double(imread([InputDir FileNames{nn}],'PixelRegion',{[RowFirst RowLast], [1 dimy]}))',HotPixThresh);
    end
end
% Print info
fprintf('Read rows %u to %u of %u %s-images of dimension [%u x %u] into sinogram of dimension [%u x %u x %u] in %fs ',...
    RowFirst,RowLast,NumImRead,ImFormat,dimx,dimy,size(sino,1),size(sino,2),size(sino,3),toc)
% Save sinogram
tic
if ~isempty(OutputDir)
    CheckTrailingSlash(OutputDir);
    CheckAndMakePath(OutputDir);
    switch OutputFormat
        case 'tif'
            for nn = 1:numel(RowsToRead)
                write32bitTIF(sprintf('%ssino_%04u.tif',OutputDir,RowsToRead(nn)),sino(:,:,nn));
            end
        case 'edf'
            for nn = 1:numel(RowsToRead)
                edfwrite(sprintf('%ssino_%04u.edf',OutputDir,RowsToRead(nn)),sino(:,:,nn)','float32');
            end
    end
    fprintf('and wrote sinos ')
end
fprintf('in %fs.\n',toc)