function im = norm_mean_std(im)
% Normalize 2D image by mean subtraction and division by the standar
% deviation.
% ARGUMENTS:
%   im: 2D matrix
% RETURNS:
%   im: 2D matrix

im = (im  - mean2(im)) / std2(im);