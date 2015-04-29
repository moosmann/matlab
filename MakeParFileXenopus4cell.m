function MakeParFileXenopus4cell(StartEndVoxels,NumOfFirstAndLastProjection,IMAGE_PIXEL_SIZE,ANGLE_BETWEEN_PROJECTIONS,ParentPath)
% For details see function 'MakeParFile'.

%% Default arguments
if nargin < 1
    StartEndVoxels{1} = [ 263  61 31;1838 1794 1524];
    StartEndVoxels{2} = [221 17 31; 1798 1746 1507];
end
if nargin < 2
    NumOfFirstAndLastProjection = [1 1600 1599];
end
if nargin < 3
    IMAGE_PIXEL_SIZE = 0.745;
end
if nargin < 4
    ANGLE_BETWEEN_PROJECTIONS = 360/1599;
end
if nargin < 5
    %ParentPath = '/mnt/tomoraid3/user/moosmann/Xenopus_4cell/';
    ParentPath = '/mnt/tomoraid-LSDF/tomo/ESRF_MI1079_ID19_July2011_inlineTomo/Xenopus_4cell';
end

MakeParFile(StartEndVoxels,NumOfFirstAndLastProjection,IMAGE_PIXEL_SIZE,ANGLE_BETWEEN_PROJECTIONS,ParentPath)