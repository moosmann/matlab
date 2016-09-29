function im = Binning(im)
% 2 x 2 binning of 2-D image. Binned pixels are normalized by factor 4.
% Images is cropped before binning such that mod(size(im),4) = 0.
% 
% Written by Julian Moosmann, last version 2013-11-08
%
% im = Binning(im)

%% MAIN
% crop to even number of pixels
im = im(1:size(im,1)-mod(size(im,1),4),1:size(im,2)-mod(size(im,2),4));
% bin
im = (im(1:2:end,1:2:end)+im(2:2:end,1:2:end)+im(1:2:end,2:2:end)+im(2:2:end,2:2:end))/4;
