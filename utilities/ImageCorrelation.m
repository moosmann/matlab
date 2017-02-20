function [out, CorMap] = ImageCorrelation(im1, im2, printInfo, showPlots, use_gpu)
%Find the relative movement of image im2 w.r.t. to im1 by means of
%correlating the images. Thus rotation axis can be determined. Allows for
%subpixel precsision.
%
% Written Julian Moosamnn. Last modification: 2017-02-20, GPU support added

% TODO: fix plotting option if axis is close to boundary and center of mass
% region is too large

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 3
    printInfo = 0;
end
if nargin < 4
    showPlots = 0;
end
if nargin < 5
    use_gpu = 0;
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Halfwidth of the region around the maximum of the correlation map which
% is used to compute the center of mass (COM) around this maximum peak.
COMwidth = 4;
[dimx, dimy] = size(im1);
% Center of images
CenX = (dimx+1)/2;
CenY = (dimy+1)/2;

% Compute the correlation map of the two images.
if ~use_gpu
    CorMap = fftshift(ifft2(fft2(im1).*fft2(rot90(im2,2))));
else    
    CorMap = fftshift( ifft2( fft2(gpuArray(im1)).*fft2(rot90(gpuArray(im2),2)), 'symmetric') );
    clear im1 im2;
end

% Find the value and the (index) position of the maximum of the correlation map.
[CorMapMaxVal, CorMapMaxInd] = max(CorMap(:));

% Convert linear index to row and column subscripts
[CorMaxPosX, CorMaxPosY] = ind2sub([dimx dimy],CorMapMaxInd);

% Add 0.5 or 1 pixel to peak position for even or odd dimensions
switch mod(dimx,2)
    case 0
        offsetX = 0.5;
    case 1
        offsetX = 1;
end
switch mod(dimy,2)
    case 0
        offsetY = 0.5;
    case 1
        offsetY = 1;
end

% Check if border of region which is used to compute COM is larger than image
if CorMaxPosX > COMwidth && CorMaxPosY > COMwidth && (dimx-CorMaxPosX) > COMwidth && (dimy-CorMaxPosY) > COMwidth
    % Region around the peak of the correlation map used to compute the COM of
    % this peak.
    CorMapRegion    = CorMap(CorMaxPosX+(-COMwidth:COMwidth),CorMaxPosY+(-COMwidth:COMwidth));
    % Position of the COM around the peak of the correlation map.
    % x-,y-coordinate-matrices needed to compute COM of the correlation map.
    [mdy, mdx]   = meshgrid(1:dimy,1:dimx);
    CorMapRegionSum = sum(CorMapRegion(:));
    comX = sum(sum(mdx(CorMaxPosX+(-COMwidth:COMwidth),CorMaxPosY+(-COMwidth:COMwidth)).*CorMapRegion))/CorMapRegionSum+offsetX;
    comY = sum(sum(mdy(CorMaxPosX+(-COMwidth:COMwidth),CorMaxPosY+(-COMwidth:COMwidth)).*CorMapRegion))/CorMapRegionSum+offsetY;
else
    comX = CorMaxPosX + offsetX;
    comY = CorMaxPosY + offsetY;
end
CorMaxPosX = CorMaxPosX + offsetX;
CorMaxPosY = CorMaxPosY + offsetY;
% Pixel shift of image 'im2' w.r.t. image 'im1'.
out.Xshift = gather( comX - CenX );
out.Yshift = gather( comY - CenY );
out.HorizontalRotationAxisPosition = gather( CenX/2 + comX/2 );
out.VerticalRotationAxisPosition   = gather( CenY/2 + comY/2 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Print results
if printInfo
    fprintf('Image dimensions: %u x %u. (MATLAB notation: dimX x dimY, x-axis directing downwards.)\n',dimx,dimy);
    fprintf('Image center: [%g %g].\n',CenX,CenY);
    fprintf('Correlation map peak: [%.2f %.2f], value: %g.\n',CorMaxPosX,CorMaxPosY,CorMapMaxVal);
    fprintf('Center of mass in the region +/- %u pixels around peak: [%f %f].\n',COMwidth,comX,comY);
    fprintf('Vertical and horizontal shift of 2nd input image w.r.t. 1st input image: [%f %f].\n',out.Xshift,out.Yshift);
    fprintf('Rotation axis position: %f (vertical rotation), %f (horizontal rotation)\n',out.VerticalRotationAxisPosition,out.HorizontalRotationAxisPosition);
end
%% Show plots
if showPlots
    % Halfwidth of the region around the computed correlation maximum which is
    % plotted.
    plotwidth = 14;
    CorMap = gather( abs( CorMap ) );
    imtool(CorMap,[],'InitialMagnification','fit')
    imtool(CorMap(CorMaxPosX+(-plotwidth:plotwidth),CorMaxPosY+(-plotwidth:plotwidth)),[],'InitialMagnification','fit')
end
