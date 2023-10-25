function im = norm_mean_std_max(im)
% Normalize 2D image by mean subtraction and division by the standar
% deviation and subsequent division by the maximum value.
% ARGUMENTS:
%   im: 2D matrix
% RETURNS:
%   im: 2D matrix

im = (im  - mean2(im)) / std2(im);
im = im / max2(im);