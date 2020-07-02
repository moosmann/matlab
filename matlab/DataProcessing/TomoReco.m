function im = TomoReco(sino,RotAxisPos,OverallAngle,OutputSize,im1,im2)

%% Default arguments.
if nargin < 2
    RotAxisPos = 528.5;
end
if nargin < 3
    OverallAngle = 180;
end
if nargin < 4
    % 0 (default): dimx/sqrt(2), 1: dimx, >1: OutputSize, <0: abs()*dimx
    OutputSize = -0.5;
end
%% Function body.
% if more than images than the need projections for tomographic
% reconstruction were taken
ProjDecr = 5;
[dim1 dim2] = size(sino);
sinoCen = dim1/2;
% Number of projections.
NumProj = dim2 - ProjDecr;
% Angle increment.
AngleIncr = OverallAngle/NumProj;
% Position of axis of rotation
if RotAxisPos == 0
    out = ImageCorrelation(im1,im2,0,0);
    RotAxisPos  = out.VerticalRotationAxisPosition;
end
% Crop sino that its center coincides with the rotation axis
% ESRF data: images to correlate are NumProj+1 and NumProj+3
pixshift = abs(round(2*(sinoCen-RotAxisPos)));
if round(RotAxisPos) < round(sinoCen)
    sino     = sino(1:end-pixshift,:);
elseif round(RotAxisPos) > round(sinoCen)
    sino     = sino(1+pixshift:end,:);
end
[dim1crop dim2crop] = size(sino);
%% Tomographic reconstruction.
% INTERPOLATION METHOD:
% 'nearest': Nearest-neighbor interpolation
% 'linear': Linear interpolation (the default)
% 'spline': Spline interpolation
% 'pchip': Shape-preserving piecewise cubic interpolation
% 'cubic': Same as 'pchip'
% FILTER to use for frequency domain filtering:
% 'Ram-Lak': Cropped Ram-Lak or ramp filter. This is the default. The frequency response of this filter is | f |. Because this filter is sensitive to noise in the projections, one of the filters listed below might be preferable. These filters multiply the Ram-Lak filter by a window that deemphasizes high frequencies.
% 'Shepp-Logan': Multiplies the Ram-Lak filter by a sinc function
% 'Cosine': Multiplies the Ram-Lak filter by a cosine function
% 'Hamming': Multiplies the Ram-Lak filter by a Hamming window
% 'Hann': Multiplies the Ram-Lak filter by a Hann window
% 'None': No filtering. When you specify this value, iradon returns unfiltered backprojection data.
% OUTPUT SIZE:
% output_size is a scalar that specifies the number of rows and columns in
% the reconstructed image. If output_size is not specified, the size is
% determined from the length of the projections: output_size = 2*floor(size(R,1)/(2*sqrt(2)))
output_size_default = 2*floor(dim1crop/(2*sqrt(2)));
if OutputSize == 0
    % MATLAB'S default
    output_size = output_size_default;
elseif OutputSize == 1;
    output_size = dim1crop;
elseif OutputSize < 0
    output_size = abs(OutputSize)*dim1crop;
else
    output_size = OutputSize;
end
output_size = round(output_size);
disp(output_size)
% Inverse Radon Transformation
im = iradon(sino,AngleIncr,'nearest','Ram-Lak',1,output_size);
fprintf('Dimension of input sinogram: %u x %u\n',dim1,dim2)
fprintf('Number of projections: %u\n',NumProj)
fprintf('Rotation axis position: %.1f\n',RotAxisPos)
fprintf('Center of sinogram: %.1f\n',sinoCen)
fprintf('pixel shift: %.1f\n',pixshift)
fprintf('Dimension of sinogram after cropping: %u x %u\n',dim1crop,dim2crop)
fprintf('Default output size: %u\n',output_size_default)
fprintf('Dimensions of reconstructed slice: %u x %u\n',size(im))
%% Show images.
%itool(im);
itool(RemoveLowFreq(im))
