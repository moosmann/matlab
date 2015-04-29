function CheckShift180(FilePrefix,NumberOfProjections,CropArea)
% Load flat-and-dark-field corrected images (int_*.edf) to check for
% vertical shift (and determine rotation axis). Work only for ESRF like
% file structure with additional quali images taken after the scan.

%% Defaul arguments
if nargin < 1
    FilePrefix = 'int';
end
if nargin < 2
    NumberOfProjections = 1599;
end
if nargin < 3
    CropArea = 0;
end

% First projection: 0 degree
imFirst = pmedfread(sprintf('%s_0001.edf',FilePrefix))';
[dimx dimy] = size(imFirst);
fprintf('  Image dimensions: [DimX DimY] = [%u %u] (MATLAB notation: X = vertical, Y = horizontal\n',dimx,dimy);
% 90 degree
imMedial = pmedfread(sprintf('%s_%04u.edf',FilePrefix,round(NumberOfProjections/2)))';
% last projection: 180 degree - 1*angle_increment
imLast = fliplr(pmedfread(sprintf('%s_%04u.edf',FilePrefix,NumberOfProjections))');
% 180 degree
imLastPlus1 = fliplr(pmedfread(sprintf('%s_%04u.edf',FilePrefix,NumberOfProjections+1))');
% 90 degree
imLastPlus2 = pmedfread(sprintf('%s_%04u.edf',FilePrefix,NumberOfProjections+2))';
% 0 degree
imLastPlus3 = pmedfread(sprintf('%s_%04u.edf',FilePrefix,NumberOfProjections+3))';
if CropArea(1) ~= 0
    x = CropArea(1,1):CropArea(1,2);
    y = CropArea(2,1):CropArea(2,2);
    yflipped = fliplr(dimy-(y));
    imFirst     = imFirst(x,y);
    imMedial    = imMedial(x,y);
    imLast      = imLast(x,yflipped);
    imLastPlus1 = imLastPlus1(x,yflipped);
    imLastPlus2 = imLastPlus2(x,y);
    imLastPlus3 = imLastPlus3(x,y);
    fprintf('  Cropped Area: y = [%u %u]; yflipped = [%u %u]\n',y(1),y(end),yflipped(1),yflipped(end));
end
ishow([imFirst imMedial imLast;imLastPlus1 imLastPlus2 imLastPlus3])
% Following three correlation should be very similiar in horizontal shift
% (rotation axis)
fprintf(1,'\nHORIZONTAL SHIFT DETERMINES ROTATION AXIS:\n');
fprintf(1,'\n1st projection at 0 degree and 1st quali image at 180 degree:\n');
ImageCorrelation(imFirst,imLastPlus1);
fprintf(1,'\n1st quali image at 180 degree and 3rd quali image at 0 degree:\n');
ImageCorrelation(imLastPlus3,imLastPlus1);
fprintf(1,'\n1st projection at 0 degree and last projection at 180 degree - 1*angle_increment:\n');
ImageCorrelation(imFirst,imLast);
% Vertical Shift can be determined. Horizontal should be zero since both
% images it's were taken at the starting angle (0 degree)
fprintf(1,'\nDETERMINE VERTICAL SHIFT:\n');
fprintf(1,'\nFirst projection at 0 degree and third quali image at 0 degree:\n');
ImageCorrelation(imFirst,imLastPlus3);
fprintf(1,'\nMedial projection at 90 degree and second quali image at 90 degree:\n');
ImageCorrelation(imMedial,imLastPlus2);

