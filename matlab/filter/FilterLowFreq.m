function [im, filtStr] = FilterLowFreq(im,sigma,hsize,BoundaryOption)
% Substract large scale structure (low frequencies) from image by means of
% blurring the image with an Gaussian low pass filter and then subtracting
% that image.
%
% sigma: scalar specifying (positive) standard deviation of the Gaussian
%        low pass filter of size 'hsize', if 0: return input image 'im' and
%        empty string 'filtStr'
% hsize: vector specifying the number of rows and columns in the Gaussian
%        low pass filter.
% BoundaryOption: Input image values outside the bounds of the image.
%
% Written by Julian Moosmann, last version: 2013-10-24
%
% im = FilterLowFreq(im,sigma,hsize,BoundaryOption)
    
    
%% Defaults arguments
if nargin < 2
   sigma = 40;
end
if nargin < 3
   hsize = sigma*[3 3];
end
if nargin<4
    BoundaryOption = 'replicate';
end
%% 
if sigma == 0
    filtStr = '';
    return;
end

%% MAIN
filtStr = sprintf('bgsS%03uV%03uH%03u',sigma,hsize);
im = im-imfilter(im,fspecial('gaussian',hsize,sigma),BoundaryOption);

