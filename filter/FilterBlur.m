function im = FilterBlur(im,SizeOfMedfiltMask_XY,SizeOfDisk)
% Blur image.
%
% Written by Julian Moosmann, 2013-10-21

%% Defaults %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    SizeOfMedfiltMask_XY = [3 3];
end
if nargin < 3
    SizeOfDisk = 10;
end
%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if SizeOfMedfiltMask_XY(1) > 0
    im = medfilt2(im,SizeOfMedfiltMask_XY,'symmetric');
end
if SizeOfDisk > 0
    %mask = imfilter(mask,fspecial('gaussian',[3 3],10),'symmetric');
    im = imfilter(im,fspecial('disk',SizeOfDisk),'symmetric');
end
