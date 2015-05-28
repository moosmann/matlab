function [sino,sinoFiltStr] = FilterSino(sino,NumAdjLines,OrientationOfLine,mask)
% Filtering of input sinogram 'sino' in order to remove ring artifacts in
% the subsequent tomographic reconstruction. Take care to choose the proper
% direction 'OrientationOfLine' in Fourier space. Output: Filtered
% sinogram. Optional output: sring containing the parameters used.
%
% NumAdjLines : scalar. Thickness of the band which is used for median
% computiation and which is centred about the line to be filtered (in
% Fourier space).
%
% OrientationOfLine : scalar/char. Default: 1.
%   Variants:
%    0: no filtering is done. 
%    1: vertical Matlab direction (1st dimension)
%    2: horizontal Matlab direction (2nd Matlab dimension)
%    3: both directions
%    'triangle': see optional argument 'mask'
%
% mask : matrix of dimensions size(sino). Optional. Default: None. Instead
% of median filter use filter mask created by MaskingTriangle.
%
% (Function formerly called 'FiltSino'.)
%
% Written by Julian Moosmann, 2013-10-18. Last modification: 2015-05-28
%
% [sino,sinoFiltStr] = FilterSino(sino,NumAdjLines,OrientationOfLine,mask)

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    NumAdjLines = 2;
end
if nargin < 3
    OrientationOfLine(1) = 1;
end
%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if OrientationOfLine == 0
    sinoFiltStr = '';
    return
end
% FFT
sino = fft2(sino);
% replace line at zero frquency by the median of the band containing this
% line at the centre
switch lower(OrientationOfLine)
    case {1,'vertical'}
        sino(:,1) = median(sino(:,[1:2+NumAdjLines end-NumAdjLines:end]),2);
        sinoFiltStr = '_filtSinoVert';
    case {2,'horizontal'}
        sino(1,:) = median(sino([1:2+NumAdjLines end-NumAdjLines:end],:),1);
        sinoFiltStr = '_filtSinoHorz';
    case {3,'both'}
        sino(:,1) = median(sino(:,[1:2+NumAdjLines end-NumAdjLines:end]),2);
        sino(1,:) = median(sino([1:2+NumAdjLines end-NumAdjLines:end],:),1);
        sinoFiltStr = '_filtSinoBoth';
    case {4,'triangle'}
        sino = mask.*sino;
        sinoFiltStr = '_filtSinoTriangle';
end
% iFFT
sino = real(ifft2(sino));