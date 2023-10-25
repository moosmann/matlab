function im = norm_max(im)
% Normalize 2D image by division of the maximum.
% ARGUMENTS:
%   im: 2D matrix
% RETURNS:
%   im: 2D matrix

im = im / max2(im);