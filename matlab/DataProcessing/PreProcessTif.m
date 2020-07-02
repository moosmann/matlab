function PreProcessTif(MedFilterThreshold,PixelRegion,NumDarks,NumFlats,NumRadios,tiffile)
 
if nargin < 1
    MedFilterThreshold = 1.1;
end
if nargin < 2
    PixelRegion =   [351 1872; 101 3908];
end
if nargin < 3
    NumDarks = 5;
end
if nargin < 4
    NumFlats = 3;
end
if nargin < 5
    NumRadios = 1;
end
if nargin < 6
    tiffile = 'radio.tif';
end

% PROGRAMME.    
% Read multi tif file. Read
stack = readmultitif(tiffile)/4.0;
if PixelRegion > 0
    stack = stack(PixelRegion(1,1):PixelRegion(1,2),PixelRegion(2,1):PixelRegion(2,2),:);
    fprintf('Size of cropped stack: %u %u %u\n',size(stack))
end
%assignin('base',['stack'],stack);
[dimx dimy NumTotal] = size(stack);
NumPos  = (NumTotal-NumDarks)/(1+NumFlats);
% Median of Dark fields. Apply hot-pixel filter to median filtered dark.
dark   = filter_im(median(stack(:,:,1:NumDarks),3),MedFilterThreshold);
int = zeros(dimx,dimy,NumPos);
for n=1:NumPos,
    int(:,:,1) = (stack(:,:,NumDarks+(NumFlats+NumRadios)*n) - dark)./ ... 
        (filter_im(median(stack(:,:,NumDarks+(NumFlats+NumRadios)*(n-1)+(1:NumFlats)),3),MedFilterThreshold) - dark);
end;

assignin('base','int',int);
