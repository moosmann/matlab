function [sino,sinoFiltStr] = FilterSino(sino,NumAdjLines,OrientationOfLine,mask)
% Filtering of input sinogram 'sino'. 
%
% NumAdjLines: scalar. Thickness of band which is used for median
% computiation and which centrally contains the line to be filtered (in FS)
%
% OrientationOfLine: scalar/char. 0, no filtering is done. 1 (default): vertial Matlab
% direction (1st dimension), 2: horizontal Matlab direction (2nd Matlab
% dimension), 3: both directions
% mask: to be created by MaskingTriangle
% Formerly known as 'FiltSino' due to existing script of the same name.
%
% Written by Julian Moosmann, 2013-10-18

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