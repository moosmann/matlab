function Xshift = LineCorrelation(line1,line2,printInfo,showPlots)
%Find the relative movement of line2 w.r.t. to line1 by means of
%correlating the lines.

%% Default arguments.
if nargin < 3
    printInfo = 1;
end
if nargin < 4
    showPlots = 0;
end

% Convert to column vectors: [dimx x 1]
line1 = line1(:);
line2 = line2(:);
%% Parameters.
% Halfwidth of the region around the maximum of the correlation map which
% is used to compute the center of mass (COM) of this maximum peak.
COMwidth = 4;
% Halfwidth of the region around the computed correlation maximum which is
% plotted.
plotwidth = 14;
%% Programme.
[dimx dimy] = size(line1);
% x-,y-coordinate-matrices need to compute COM of the correlation map.
[~, mdx]   = meshgrid(1:dimy,1:dimx);
% Compute the correlation map of the two lines.
CorMap      = fftshift(ifft(fft(line1).*fft(rot90(line2,2))));
% Find the value and the (index) position of the maximum of the correlation map.
[CorMapMaxVal CorMapMaxInd] = max(CorMap(:));
% Convert linear index to row and column subscripts
[CorMaxPosX] = ind2sub([dimx dimy],CorMapMaxInd);
% Region around the peak of the correlation map used to compute the COM of
% this peak.
CorMapRegionRange = 1+mod(-1+(CorMaxPosX+(-COMwidth:COMwidth)),dimx);
CorMapRegion    = CorMap(CorMapRegionRange);
CorMapRegionSum = sum(CorMapRegion(:));
% Position of the COM around the peak of the correlation map.
comX = sum(sum(mdx(CorMapRegionRange).*CorMapRegion))/CorMapRegionSum;
%comY = sum(sum(mdy(CorMaxPosX+(-COMwidth:COMwidth),CorMaxPosY+(-COMwidth:COMwidth)).*CorMapRegion))/CorMapRegionSum;
% Pixel shift of line 'line2' w.r.t. line 'line1'.
Xshift = comX-dimx/2;
%out.Yshift = comY-dimy/2;
%out.HorizontalRotationAxisPosition = dimx/4+comX/2;
%out.VerticalRotationAxisPosition   = dimy/4+comY/2;
%% Print results.
if printInfo
    fprintf('line dimensions: %u x %u. (MATLAB notation: dimX x dimY, x-axis directing downwards.)\n',dimx,dimy);
    fprintf('Center of lines: [%g %g].\n',dimx/2);
    fprintf('Maximum of correlation map found at [%u] with a value of %g.\n',CorMaxPosX,CorMapMaxVal);
    fprintf('Center of mass in the region +/- %u pixels around maximum: [%f].\n',COMwidth,comX);
    fprintf('Vertical and horizontal shift of 2nd input line w.r.t. 1st input line: [%f].\n',Xshift);
    %fprintf('Rotation axis position: %6f (horizontal rotation)\n',out.HorizontalRotationAxisPosition);
end
%% Print plots.
if showPlots
    %imtool(CorMap,[],'InitialMagnification','fit')
    imtool(CorMap(CorMaxPosX+(-plotwidth:plotwidth)),[],'InitialMagnification','fit')
end
