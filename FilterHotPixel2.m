function [imFiltered,imMedian,mask,R,NumOfFilteredPixels] = FilterHotPixel2(im,threshold,printInfo,MedianFilterRadius,filterDeadPixel,GaussianSigma,imMedian)
%Script is written to remove salt and pepper noise from the input data by
%the median value. By default 1% of the pixel will be filtered. 
%FilterHotPixel2 provides additonal functionalities compared to
%FilterHotPixel.
%Written by J. Moosmann, Jul 2011, modified Nov 2012, code deprecated
%compared to FilterHotPixel2
%
% im: image to be filterd
% threshold: if < 1, interpreted as percentage of pixel which should be
% filtered, default: 1%. if > 1, value for the matrix of the ratio
% im/median(im) above what pixel will be filtered. 
% MedianFilterRadius: Radii of the median filter applied to the image. This
% median filtered image is needed for the ratio matrix which determines the
% pixel to be filtered
% GaussianSigma: Gaussian filter which can be applied to the input im
% before building the ration im/median(im), thus gaussian(im)/median(im)
% imMeadian: = median(im), default: if empty will be generated with above
% given values for the radii
%
%[imFiltered,imMedian,mask,R,NumOfFilteredPixels] = filter_im(im,threshold,MedianFilterRadius,GaussianSigma,imMedian)

%% Default values.
if nargin < 2
    threshold = 0.01;% 1% of the pixel will be filtered
end
if nargin < 3
    printInfo = 0;
end
if nargin < 4
    MedianFilterRadius = [3 3];
end
if nargin < 5
    filterDeadPixel = 1;
end
if nargin < 6
    GaussianSigma = 0;
end
% Return if threshold is zero.
if threshold == 0
    imFiltered = im;
    return
end
%% Parameters.
medx = MedianFilterRadius(1,1);
medy = MedianFilterRadius(1,2);
NumOfPixels = numel(im);
%% Median filterd image. 
%If not given as input it is generated. Note the option 'symmetric' in
%medfilt2 to avoid zero values in the median filtered image.
if nargin < 7
    imMedian = medfilt2(im,[medx medy],'symmetric');
    % imMedian = medfilt2(padarray(im,[medx medy],'symmetric','both'),[medx medy]);
    % imMedian = imMedian(1+medx:end-medx,1+medy:end-medy);
end
%% Ratio of image and median filtered image
% optional Gaussian filtering of nominator before building ratio
switch GaussianSigma > 0
    case 1
        R = imfilter(im,fspecial('gaussian',ceil(2*GaussianSigma*[medx medy]),GaussianSigma))./imMedian;
    case 0
        R = im./imMedian;
end
imFiltered = im;
%% Filter dead pixel
if filterDeadPixel > 0
    mask = imFiltered == 0;
    imFiltered(mask) = imMedian(mask);
    NumDead = numel(mask(mask));
end
%% Determine threshold if not given explixitly as a value > 1.
if threshold < 1
    Rsorted = sort(R(:));
    threshold = Rsorted(floor(NumOfPixels*(1-threshold)));
end
%% Create binary matrix of pixels which are median filtered.
mask = R > threshold;
%% Get filtered image
% Assigen to the pixels of R which are above the threshold, the
% median filtered value.
imFiltered(mask) = imMedian(mask);
NumOfFilteredPixels = numel(imFiltered(mask));
%% Print info.
if printInfo
    [imMin imMax imMean imStd] = Domain(im,1,'unfiltered image',0);
    [imFilteredMin imFilteredMax imFilteredMean imFilteredStd] = Domain(imFiltered,1,'  filtered image',0);
    fprintf(1,'Pixels filtered: %u of %u (%3.2f%%). ',NumOfFilteredPixels,NumOfPixels,100*NumOfFilteredPixels/NumOfPixels);
    fprintf(1,'[Min Max Mean Std]=[%g %g %g %g](unfiltered)',imMin,imMax,imMean,imStd);
    fprintf(1,',[%g %g %g %g](filtered)\n',imFilteredMin,imFilteredMax,imFilteredMean,imFilteredStd);
    if filterDeadPixel > 0
        fprintf('Number of dead pixels in input image: %u, in filtered image %u. Minimum of filtered image %u.\n',...
            NumDead,numel(imFiltered(imFiltered==0)),imFilteredMin)
    end
end
end

    
