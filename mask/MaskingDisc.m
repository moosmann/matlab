function im = MaskingDisc(im,radFac,maskVal)
% Keep central disc with radius 'radFac'*size(im,1) and set the region
% outside of the disk to the mean within the disc.
%
% im: 2D-matrix
% radFac: scalar in [0 1]. default 0.95. defines the radius of the disc
%
% Written by Julian Moosmann, last version2014-01-30
%
% im = MaskingDisc(im,radFac) 

if nargin < 2
    radFac = 0.95;
end
if nargin < 3
    maskVal = [];
end

[x,y] = meshgrid(-size(im,2)/2:size(im,2)/2-1,-size(im,1)/2:size(im,1)/2-1);

m = sqrt(x.^2+y.^2) < radFac*size(im,1)/2;

if isempty(maskVal)
    % set maskVal to inner mean
    maskVal = sum(im(:).*m(:))/sum(m(:));
end

im = m.*im + (1-m)*maskVal;