function imstack = readMultiTif(filename,PixelRegion,PreFactor)

% imstack = readmultitif(filename,PixelRegion)
% Fast importing of a stack of images form a multitif file.
%
% PixelRegion is given by cell struct: {[dim_x_first dim_x_last]
% [dim_y_first dim_y_last]}, where dim_x is the vertical 
% and dim_y the horizontal direction for arrays in MATLAB. Default is
% whole image.

%  Preallocation of memory speeds up importing tremedously.

% Get information about each image in the stack included in the graphics file.
tinfo = imfinfo(filename,'tif');
% Get number of images.
NumOfImages = numel(tinfo);
% Get dimensions of first image. (Assuming that all images have the same dimensions.
dimx = tinfo(1).Height;
dimy = tinfo(1).Width;
% Preallocate memory.
imstack = cell(1,NumOfImages);
% Import image stack.
if nargin<2
    PixelRegion={[1 dimx] [1 dimy]};
end
if nargin<3
    PreFactor = 1;
end
tic
for k = 1:NumOfImages
    imstack{k} = PreFactor*double(imread(filename,k,'PixelRegion',PixelRegion));
end
tread=toc;
% Print dimensions of multitif stack.
fprintf(sprintf(['Read multitif file of format: %i x %i x %i. Elapsed ' ...
                 'time: %3.3fs.\n'],dimx,dimy,NumOfImages,tread));
