function imStack = Readstack_GI(InputPath,FilePrefix,DataType,FirstFileToRead,StepSize,NumSteps,HotPixelFilter)
% Read multiple images into a 3d stacks. Arguments to define the images to
% read are suited for GI experiments.
%
% InputPath: String. Path to the tif files.
% FilePrefix: String. Default ''. Optional.
% DataType: String. Default: 'edf'. Image format: 'edf', 'tif', or any
% other format MATLAB's imread can handle.
% FirstFileToRead: Integer. Default 1. Number of the first image to be read
% StepSize: Integer. Default 1. (Number of images per position of grating)
% NumSteps: Integer. Default 0: Read all images.
% HotPixelFilter: Integer Vec 1x1 or 2x1. Default 0. For < 1: percentage of pixels to be
% filered. For > 1: threshold for filering using ratio im/medfilt(im)
%
% Written by J.Moosmann and V.Altapova, first version 2010-08-20, last
% version 2013-10-16 
%
%imStack = readstack(InputPath,FilePrefix,DataType,FirstFileToRead,StepSize,NumSteps,HotPixelFilter)

%% Default arguments.
if nargin<1
    InputPath = '.';
end;
if nargin<2
    FilePrefix = '';
end;
if nargin<3
    DataType = 'edf';
end;
if nargin<4
    FirstFileToRead = 1;
end;
if nargin<5
    StepSize = 1;
end;
if nargin<6
    NumSteps = 0;
end;
if nargin<7
    HotPixelFilter = 0;
end;
%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CheckTrailingSlash(InputPath)
%% Get struct containing the file names
imStruct = dir([InputPath FilePrefix '*.' DataType]);
%% Get and set stack parameter
switch lower(DataType)
    case {'edf'}
        imStack = pmedfread([InputPath imStruct(1).name])';
    case {'tif'}
        imStack = imread([InputPath imStruct(1).name]);
end
NumIm = numel(imStruct);
[dim1,dim2]  = size(imStack);
% Define which images to read.
if NumSteps == 0,
    stepList = FirstFileToRead:StepSize:NumIm;
else
    stepList = FirstFileToRead:StepSize:(FirstFileToRead+(NumSteps-1)*StepSize);
end
NumRead = length(stepList);
%% Preallocation
imStack = zeros(dim1,dim2,NumRead);
%% Loop over images
imCounter = 1;
tic;
switch lower(DataType)
    case {'edf'}
        for kk = stepList
            imStack(:,:,imCounter) = FilterPixel(pmedfread([InputPath imStruct(kk).name])',HotPixelFilter);
            imCounter = imCounter + 1;
        end
    case {'tif'}
        for kk = stepList,
            imStack(:,:,imCounter) = FilterPixel(imread([InputPath imStruct(kk).name]),HotPixelFilter);
            imCounter = imCounter + 1;
        end
end
%% Print information 
fprintf(['Read %i of %i ''' DataType ''' images in %.2fs contained in\n' InputPath '\n'],NumRead,NumIm,toc);
if HotPixelFilter > 0
    fprintf('Hot-pixel filter applied. Filtered %g%% of all pixels.\n',HotPixelFilter*100)
end
