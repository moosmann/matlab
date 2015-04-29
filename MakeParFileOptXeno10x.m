function MakeParFileOptXeno10x(StartEndVoxels,NumOfFirstAndLastProjection,IMAGE_PIXEL_SIZE,ANGLE_BETWEEN_PROJECTIONS,ParentPath)
% For details see function 'MakeParFile'.

%% Default arguments
if nargin < 1
    %StartEndVoxels = 0;
    StartEndVoxels{1} = [ 43 83 21; 936 906 740];
    StartEndVoxels{2} = StartEndVoxels{1};
    StartEndVoxels{3} = [ 29 73 21; 932 896 740];
    StartEndVoxels{4} = StartEndVoxels{3};
end
if nargin < 2
    NumOfFirstAndLastProjection = [1 749 749];
end
if nargin < 3
    IMAGE_PIXEL_SIZE = 1.5;
end
if nargin < 4
    ANGLE_BETWEEN_PROJECTIONS = 360/749;
end
if nargin < 5
    %ParentPath = '/mnt/tomoraid3/user/moosmann/Xenopus_4cell/';
    ParentPath = '/mntdirect/_data_visitor/mi1057/bm05/pc/opt_xeno_stage12_tomo_10x/';
end

MakeParFile(StartEndVoxels,NumOfFirstAndLastProjection,IMAGE_PIXEL_SIZE,ANGLE_BETWEEN_PROJECTIONS,ParentPath)