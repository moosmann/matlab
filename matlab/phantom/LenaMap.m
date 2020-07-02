function out = LenaMap()
%% Read lena test image and normalized it.
out = zeros(1024);
out(257:768,257:768) = normat(double(imread('~/lena.tif')));
