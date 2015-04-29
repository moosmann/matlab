function im = mmm(im)
% Subtract mean from input matrix
%
%Written by Julian Moosmann, 2013-10-09

im = im - mean(im(:));