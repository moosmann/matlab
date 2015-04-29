function MakeParFileZebrafish(StartEndVoxels,NumOfFirstAndLastProjection,IMAGE_PIXEL_SIZE,ParentPath,DataPrefix)

%% Default arguments
if nargin < 1
    %StartEndVoxels{1} = [];
    StartEndVoxels{1} = [ 121 441 474;1340 1910 1030];
    StartEndVoxels{5} = [231 447 1; 1000 1110 2048];
end
if nargin < 2
    NumOfFirstAndLastProjection(1:5)  = {[1 1599 1599]};
    NumOfFirstAndLastProjection{5}    = [1 999 999];
end
if nargin < 3
    IMAGE_PIXEL_SIZE = 0.745;
end
if nargin < 4
    ANGLE_BETWEEN_PROJECTIONS(1:5) = {360/1599};
    ANGLE_BETWEEN_PROJECTIONS{5}   = 360/999;
end
if nargin < 5
    %ParentPath = '/mnt/tomoraid-LSDF/tomo/ESRF_MI1079_ID19_July2011_inlineTomo_Zebrafish';
    ParentPath = '/mnt/tomoraid-LSDF/tomo/ESRF_MI1079_ID19_July2011_inlineTomo/ZebraFish';
end
if nargin < 6 
    DataPrefix = 'ZebraFish';
end

%% Start MakeParFile loop.
MakeParFile(StartEndVoxels,NumOfFirstAndLastProjection,IMAGE_PIXEL_SIZE,ANGLE_BETWEEN_PROJECTIONS,ParentPath,DataPrefix)
