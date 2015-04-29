function [w, wx, wy, w1x, w1y] = FilterGaussian(InverseGaussianWidth,DimensionsXY)
%
%
% Written by Julian Moosmann, when ?, what does it do ??

dimx = DimensionsXY(1);
dimy = DimensionsXY(2);

%% Create Fourier space filter
w1x = 1- window(@gausswin,dimx,InverseGaussianWidth);
w1y = 1- window(@gausswin,dimy,InverseGaussianWidth);
[wx, wy] = meshgrid(w1y,w1x);
w = fftshift(wx.*wy);
w(1,1) = 1;% Restore mean