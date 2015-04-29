function MakeParFileXenopusESRF(StartEndVoxels,NumOfFirstAndLastProjection,IMAGE_PIXEL_SIZE,ParentPath)

%% Default arguments
if nargin < 1
    StartEndVoxels{1} = [303 489 197; 1236 1316 1052];
    StartEndVoxels{2} = [145 259 129; 1062 1070 1038];
    StartEndVoxels{3} = [215 191 339;  970  952 1666];
    StartEndVoxels{4} = [303 241 146;  982  922 2000];
    StartEndVoxels{5} = [265 297 143;  896 1050 2048];
end
if nargin < 2
    % [FirstProjectionToUse LastProjectionToUse NumberOfProjectionsOverFullAngle]
    % Due to PyHST problems sometimes you have to adjust the number of
    % projections used for PyHST although the number of projection used for
    % a full tomograph
    NumOfFirstAndLastProjection = [1 1600 1599];
end
if nargin < 3
    IMAGE_PIXEL_SIZE = 1.4;
end
if nargin < 4
    ANGLE_BETWEEN_PROJECTIONS = 180/1599;
end
if nargin < 5
    %ParentPath = '/mnt/tomoraid3/user/moosmann/Xenopus_ESRF_May2011/';
    %ParentPath = '/mnt/tomoraid-LSDF/tomo/Xenopus_ESRF_May2011/';
    ParentPath = '/mnt/tomoraid-LSDF/tomo/ESRF_May2011_Xenopus';
end

MakeParFile(StartEndVoxels,NumOfFirstAndLastProjection,IMAGE_PIXEL_SIZE,ANGLE_BETWEEN_PROJECTIONS,ParentPath)