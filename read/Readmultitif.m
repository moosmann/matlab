function imstack = readmultitif(filename,ImagesToRead,PixelRegion)

% imstack = readmultitif(filename,PixelRegion)
% Fast importing of a stack of images form a multitif file.
%
% PixelRegion is given by cell struct: {[dim_x_first dim_x_last]
% [dim_y_first dim_y_last]}, where dim_x is the vertical 
% and dim_y the horizontal direction for arrays in MATLAB. Default is
% whole image.
%
% Preallocation of memory speeds up importing tremedously.

if nargin < 1
    filename = 'radio.tif';
end
if nargin < 2
    ImagesToRead = 0;
end

% Get information about each image in the stack included in the graphics file.
tinfo = imfinfo(filename,'tif');
% Get number of images.
if ImagesToRead == 0
    ImagesToRead = 1:numel(tinfo);
end
NumOfImages = numel(ImagesToRead);
% Get dimensions of first image. (Assuming that all images have the same dimensions.
dimx = tinfo(1).Height;
dimy = tinfo(1).Width;
% Preallocate memory and import image stack.
warning('off','MATLAB:tifftagsread:tiffTag:wrongTagDataFormat')
if nargin < 3
    imstack = zeros(dimx,dimy,NumOfImages);
    tic;
    for k = 1:NumOfImages
        imstack(:,:,k) = double(imread(filename,ImagesToRead(k)));
    end;
    tread=toc;
else
    imstack = zeros(1+PixelRegion{1}(2)-PixelRegion{1}(1),1+PixelRegion{2}(2)-PixelRegion{2}(1),NumOfImages);
    tic;
    for k = 1:NumOfImages
        imstack(:,:,k) = double(imread(filename,ImagesToRead(k),'PixelRegion',PixelRegion));
    end;
    tread=toc;
end
% Print dimensions of multitif stack.
fprintf(sprintf(['Read multitif file. Format: %i x %i x %i. Elapsed ' ...
                 'time: %3.3fs.\n'],dimx,dimy,NumOfImages,tread));
