function im = SubtractMean(im)
% Output is the mean subtracted input.

im = im - mean(im(:));