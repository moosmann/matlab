function im = RotAxisSymmetricPadding(im,RotAxisPos,printInfo)
% Pad and/or crop image such that the rotation axis is centred at dim2/2.
%
% Written by Julian Moosmann, last modified: 2013-09-09
% im = RotAxisSymmetricPadding(im,RotAxisPos,printInfo)

%% Default arguments
if nargin < 3
    printInfo = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if RotAxisPos == 0
    return;
end
%% Parameters
RotAxisPos = round(double(RotAxisPos));
dimHor     = size(im,2);
imCenter   = dimHor/2;
OffCenter  = abs(imCenter-RotAxisPos);
ClosestPowerOf2 = 2^nextpow2(imCenter - OffCenter);
OffAxis = abs(ClosestPowerOf2-RotAxisPos);
RightResidual = abs((dimHor - RotAxisPos) - ClosestPowerOf2);
%% Padding
if RotAxisPos <= ClosestPowerOf2
    if dimHor - RotAxisPos > ClosestPowerOf2
        %fprintf('%g\n',OffAxis);
        im = padarray(im(:,1:floor(RotAxisPos+ClosestPowerOf2)),[0 OffAxis],'symmetric','pre');    
    else            
        im = padarray(im,[0 RightResidual],'symmetric','post');
        im = padarray(im,[0 OffAxis],'symmetric','pre');
    end
else
    im = padarray(im(:,1+OffAxis:end),[0 RightResidual],'symmetric','post');
end
%% Print info
if printInfo
    fprintf('RotAxisPos: %f\n',RotAxisPos);
    fprintf('Input image dimensions: [%g %g]\n',size(im));
    fprintf('Input image center: %g\n',imCenter);
    fprintf('RotAxis off center: %f\n',OffCenter);
    fprintf('Closest power of 2 of RotAxis: %f\n',ClosestPowerOf2);
    fprintf('Closest power of 2 off RotAxis: %f\n',OffAxis);
    fprintf('dimHor - RotAxisPos: %f\n',dimHor-RotAxisPos);
    fprintf('Right residual: %f\n',RightResidual);
    fprintf('Output image dimensions: [%g %g]\n',size(im));
end