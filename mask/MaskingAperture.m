function mask = MaskingAperture(dark,flat,ThreshFac)
% Create mask from dark and flat field to mask low count region due to
% aprture of optical system

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
darkMean = mean(dark(:));
% Make mask
mask = ones(size(dark));
mask(flat < ThreshFac*darkMean) = 0;
mask = medfilt2(mask,[3 3],'symmetric');
%mask = imfilter(mask,fspecial('gaussian',[3 3],10),'symmetric');
mask = imfilter(mask,fspecial('disk',10),'symmetric');

%mask = floor(mask);