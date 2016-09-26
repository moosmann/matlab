function CheckShift360(FilePrefix,NumberOfProjections,CropArea)
% Load flat-and-dark-field corrected images (int_*.edf) to check for
% vertical shift (and determine rotation axis). Work only for ESRF like
% file structure with additional quali images taken after the scan.
% 
% image		angle/??
% 1		    0
% 1599		360-1*(angle increment)
% 1600		360
% 1601		270
% 1602		180
% 1603		90
% 1604		0

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

%% Reading (and flipping if necessary.)
% first projection: 0??
imFirst = pmedfread(sprintf('%s_0001.edf',FilePrefix))';
[dimx dimy] = size(imFirst);
% medial projection: 180??
imMedial = fliplr(pmedfread(sprintf('%s_%04u.edf',FilePrefix,round(NumberOfProjections/2)))');
% last projection: 360?? - 1*angle_increment
imLast = pmedfread(sprintf('%s_%04u.edf',FilePrefix,NumberOfProjections))';
% last+1: 360??
imLastPlus1 = pmedfread(sprintf('%s_%04u.edf',FilePrefix,NumberOfProjections+1))';
% last+2: 270??
imLastPlus2 = fliplr(pmedfread(sprintf('%s_%04u.edf',FilePrefix,NumberOfProjections+2))');
% last+3: 180??
imLastPlus3 = fliplr(pmedfread(sprintf('%s_%04u.edf',FilePrefix,NumberOfProjections+3))');
% last+4: 90??
imLastPlus4 = pmedfread(sprintf('%s_%04u.edf',FilePrefix,NumberOfProjections+4))';
% last+5: 0??
imLastPlus5 = pmedfread(sprintf('%s_%04u.edf',FilePrefix,NumberOfProjections+5))';
%% Cropping.
if CropArea ~= 0
    x = CropArea(1,1):CropArea(1,2);
    y = CropArea(2,1):CropArea(2,2);
    yflipped = fliplr(dimy-(y));
    imFirst     = imFirst(x,y);
    imMedial    = imMedial(x,y);
    imLast      = imLast(x,y);
    imLastPlus1 = imLastPlus1(x,y);
    imLastPlus2 = imLastPlus2(x,yflipped);
    imLastPlus3 = imLastPlus3(x,yflipped);
    imLastPlus4 = imLastPlus4(x,y);
    imLastPlus5 = imLastPlus5(x,y);
end
%ishow([imFirst imMedial imLast;imLastPlus1 imLastPlus2 imLastPlus3])
%% Correlating.
% Following three correlation should be very similiar in horizontal shift
% (rotation axis)
fprintf(1,'\n\nHORIZONTAL AND VERTICAL SHIFT DURING SCAN:\n');
fprintf(1,'Compare images recorded at the beginning and at the end of the scan.\n');
fprintf(1,'\nProjection: 1st    at   0?? and last   at 360??-1*angle_increment:\n');
ImageCorrelation(imFirst,imLast);
fprintf(1,'\nProjection: 1st    at   0?? and last+1 at 360??:\n');
ImageCorrelation(imFirst,imLastPlus1);
fprintf(1,'\nProjection: 1st    at   0?? and last+5 at   0??:\n');
ImageCorrelation(imFirst,imLastPlus5);
fprintf(1,'\nProjection: Medial at 180?? and last+3 at 180??:\n');
ImageCorrelation(imMedial,imLastPlus3);

%% Axis of rotation.
fprintf(1,'\n\nROTATION AXIS:\n');
fprintf(1,'Compare images which are rotated about 180?? w.r.t. to each other.\n');
fprintf(1,'\nProjection: 1st    at   0?? and last+3 at 180??:\n');
ImageCorrelation(imFirst,imLastPlus3);
fprintf(1,'\nProjection: last+1 at 360?? and last+3 at 180??:\n');
ImageCorrelation(imLastPlus1,imLastPlus3);
fprintf(1,'\nProjection: last+5 at   0?? and last+3 at 180??:\n');
ImageCorrelation(imLastPlus5,imLastPlus3);
fprintf(1,'\nProjection: last+4 at  90?? and last+2 at 270??:\n');
ImageCorrelation(imLastPlus4,imLastPlus2);
