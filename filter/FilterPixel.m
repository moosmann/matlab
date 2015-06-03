function [im,HotPixPercent,DeadPixPercent,DarkPixPercent] = FilterPixel(im,FiltThreshHot_FiltThreshDark,printInfo,medianFilterRadius,filterDeadPixel,filterInfsAndNans)
%Filter hot, dark, and dead pixels.
%Substitute hot (and dead) pixels by median values. Hot and dark pixels are
%found using a ratio of the image and the median filtered images which has
%to be bigger than a FiltThreshHot_FiltThreshDark.
%
% im: matrix, image to be filtered
% FiltThreshHot_FiltThreshDark: scalar or vector, default: [0.01 0]; values < [1 0.5] are
% interpreted as the percentage of pixel to be filtered; for values > [1
% 0.5], pixels with ratio-matrix > FiltThreshHot_FiltThreshDark(1) and FiltThreshHot_FiltThreshDark(2) are filtered 
% printfInfo: scalar, default 0
% medianFilterRadius: 2-vector, default: [3 3], radii of the median filter applied to the image
% filterDeadPixel: scalar, default true
%
%If FiltThreshHot_FiltThreshDark = [0 0] and filterDeadPixel = 0 the
%input image is directly returned without doing any computations.
%
%Calculating the threshold from the prescription to filter X % of all pixel
%is more compuationally extensive compared to the case where the values are
%given directly, and thus should be avoided for a large amount of data.
%
%Printing information about the unfiltered and filtered image also produces
%unnecessary workload.
%
% Written by Julian Moosmann, first version 2011-Jul, last mod 2015-05-27
%
%[im,HotPixPercent,DeadPixPercent,DarkPixPercent] = FilterPixel(im,FiltThreshHot_FiltThreshDark,printInfo,medianFilterRadius,filterDeadPixel)


%% Default values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    FiltThreshHot_FiltThreshDark = [0.01 0];
end
if nargin < 3
    printInfo = 0;
end
if nargin < 4
    medianFilterRadius = [3 3];
end
if nargin < 5
    filterDeadPixel = 1;
end
if nargin < 6
    filterInfsAndNans = [1, 1];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check threshold argument
HotThresh  = FiltThreshHot_FiltThreshDark(1);
if length(FiltThreshHot_FiltThreshDark) == 2
    DarkThresh = FiltThreshHot_FiltThreshDark(2);
else
    DarkThresh = 0;
end
%% Return original image if Threshhold are zero.
if HotThresh == 0 && DarkThresh == 0 && filterDeadPixel == 0
    return;
end
%% Check data type for non-integerness
if isinteger(im)
    im = single(im);
end
%% Info
NumPix = numel(im);
NumDead = 0;
NumHot  = 0;
NumDark = 0;
NumNan = 0;
NumInf = 0;
if printInfo
    imMin  = min(im(:));
    imMax  = max(im(:));
    imMean = mean(im(:));
    imStd  = std(im(:));
end
mask = zeros(size(im));
% INFs
if filterInfsAndNans(1)
    mask = mask | isinf(im);
    NumInf = sum(mask(:));
end
% NANs
if filterInfsAndNans(2)
    mask = mask | isnan(im);
    NumNan = sum(mask(:)) - NumInf;
end

%% Median filterd image
% set Infs and Nans to zero before computation of median filtered map
im(mask) = 0;
imMedian = medfilt2(im, medianFilterRadius, 'symmetric');

%% Dead pixel mask
if filterDeadPixel > 0
    mask = mask | im <= 0;
    NumDead = sum(mask(:)) - NumInf -NumNan;
    if HotThresh == 0 && DarkThresh == 0 && NumDead == 0
        return
    end
end

%% Ratio of image and median filtered image
R = im./imMedian;
if (HotThresh > 0 && HotThresh <= 1) || (DarkThresh > 0 && DarkThresh <= 0.5)
    Rsorted = sort(R(:));
end
%% Hot pixel mask
if HotThresh > 0
    % Determine HotThresh if not given explixitly as a value > 1
    if HotThresh < 1
        % Set FiltThreshHot_FiltThreshDark to filter 'theshold' % of all pixels
        HotThresh = Rsorted(floor(NumPix*(1-HotThresh)));
    end
    % Add hot pixel mask
    if exist('mask','var')
        mask = mask | R > HotThresh;
    else
        mask = R > HotThresh;
    end
    NumHot = sum(mask(:)) - NumDead;
end
%% Dark pixel mask
if DarkThresh > 0
    % Determine DarkThresh if not given explixitly as a value < 1
    if DarkThresh < 0.5
        % Set FiltThreshHot_FiltThreshDark to filter 'theshold' % of all pixels
        DarkThresh = Rsorted(floor(NumPix*(DarkThresh)));
    end
    mask = mask | R < DarkThresh;
    NumDark = sum(mask(:)) - NumDead - NumHot;
end
%% Replace dead, hot, and dark pixels by median values
im(mask) = imMedian(mask);
%% Print info
if printInfo
    imFiltMin  = min(im(:));
    imFiltMax  = max(im(:));
    imFiltMean = mean(im(:));
    imFiltStd  = std(im(:));
    fprintf('Total number of pixels: %u\n',NumPix);
    if filterDeadPixel > 0
        fprintf('Number of dead pixels: %9u (%3.2f%%)\n',NumDead,100*NumDead/NumPix)
    end
    if HotThresh > 0
        fprintf('Number  of hot pixels: %9u (%3.2f%%), threshold: %9g\n',NumHot,100*NumHot/NumPix,HotThresh);
    end
    if DarkThresh > 0
        fprintf('Number of dark pixels: %9u (%3.2f%%), threshold: %9g\n',NumDark,100*NumDark/NumPix,DarkThresh);
    end
    if NumInf > 0
        fprintf('Number of INFs: %9u (%3.2f%%)\n',NumInf,100*NumInf/NumPix);
    end
    if NumNan > 0
        fprintf('Number of NANs: %9u (%3.2f%%)\n',NumNan,100*NumNan/NumPix);
    end 
    fprintf('Before: [Min Max Mean Std] = [%9g %9g %9g %9g]\n',imMin,imMax,imMean,imStd);
    fprintf('After : [Min Max Mean Std] = [%9g %9g %9g %9g]\n',imFiltMin,imFiltMax,imFiltMean,imFiltStd);
end
HotPixPercent = NumHot/NumPix;
DeadPixPercent = NumDead/NumPix;
DarkPixPercent = NumDark/NumPix;
%%% END OF CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
