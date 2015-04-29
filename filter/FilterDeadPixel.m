function [imFiltered mask] = FilterDeadPixel(im,MedianFilterSize,printInfo,imMedian)
% Remove dead pixel that are zero and replace them with the median.

%% Default arguments
if nargin < 2
    printInfo = 1;
end
if nargin < 3
    MedianFilterSize = [4 4];
end
if nargin < 4
    printInfo = 1;
end
%% Get mask of dead pixels
mask = im == 0;
% Number of dead pixels
NumDead = sum(mask(:));
%% Return if no dead pixel found
if NumDead == 0
    if printInfo
        fprintf('No dead pixels found.\n')
    end
    imFiltered = im;
    return
end
%% Median filterd image
%If not given as input it is generated. Note the option 'symmetric' in
%medfilt2 to avoid zero values in the median filtered image.
if nargin < 5
    imMedian = medfilt2(im,MedianFilterSize,'symmetric');
end
%% Replace dead pixels
imFiltered = im;
imFiltered(mask) = imMedian(mask);
%% Print info
if printInfo
    domain(im,'input image   ')
    domain(imFiltered,'filtered image')
    fprintf('Number of dead pixels in imput image   : %u\n',NumDead)
    fprintf('Number of dead pixels in filtered image: %u\n',numel(imFiltered(imFiltered==0)))
    fprintf('Minimum value of filtered image: %g\n',min(imFiltered(:)))
end
