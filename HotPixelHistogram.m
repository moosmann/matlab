function [numpix thresh] = HotPixelHistogram(im,NumBins)
% Compute a kind of histogram showing for each bin the number of pixel that
% have values above a certain values.

if nargin < 2
    NumBins = 100;
end

% Convert image to line vector to speed up computaton
im = double(im(:));
imMax = max(im);
NumPix = numel(im);
thresh = 0:imMax/(NumBins-1):imMax;

for ii = NumBins:-1:1
    m = im > thresh(ii);
    numpix(ii) = sum(m(:));
end

%figure('Name','Histogram: Absolute number of pixels above threshold VS threshold')
%plot(thresh,numpix,'.-')
figure('Name','Histogram: Percentage of pixels above threshold VS threshold')
plot(thresh,numpix/NumPix*100,'.-')

thp = input(sprintf('Enter threshold for hot pixels: '));
m = thresh > thp;
thpPos = NumBins - sum(m(:));
nhp = numpix(thpPos);
fprintf('Filtering above %g corresponds to %g filtered pixels or %f (%f%%)\n',thp,nhp,nhp/NumPix,nhp/NumPix*100)