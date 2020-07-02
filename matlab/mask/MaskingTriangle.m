function mask = MaskingTriangle(DimVert_DimHor,YoverX,offsetX_offsetY)
% Create mask of triangular shape to be multiplied in FS to filter
% frequency components in sinogram.
%
% Written by Julian Moosmann, 2013-10-18, function is a quick build,
% possibly buggy.
%
%mask = MaskingTriangle(DimVert_DimHor,YoverX,offsetX_offsetY)

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    YoverX = 210/148;
end
if nargin < 3
    offsetX_offsetY = [1 0];
end
%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creat coordinate grid
DimX = DimVert_DimHor(2);
DimY = DimVert_DimHor(1);
x = (1:DimX) - offsetX_offsetY(1);
y = (1:DimY) - offsetX_offsetY(2);
[x, y] = meshgrid(x,y);
% Create mask
mask = y-YoverX*x>0;
mask = mask & flipud(mask);
mask =  1 - (mask + fliplr(mask));
% Blur mask
%mask = imfilter(mask,fspecial('disk',5),'symmetric');
