function MakeParFileOptXeno(StartEndVoxels,NumOfFirstAndLastProjection,EffectivePixelSize,AngleBetweenProjections,ParentPath)
% For details see function 'MakeParFile'.

%% Default arguments
if nargin < 1
    StartEndVoxels = 0;
%     StartEndVoxels{1} = [ 21  103 45;1740 1760 1492];
%     StartEndVoxels{2} = StartEndVoxels{1};
%     StartEndVoxels{3} = [ 21  103 45;1740 1760 1492];
%     StartEndVoxels{4} = StartEndVoxels{3};
end
if nargin < 2
    NumOfFirstAndLastProjection(1:26)  = {[1 999]};
    NumOfFirstAndLastProjection(25:26)  = {[1 749]};
    NumOfFirstAndLastProjection(27:30)  = {[1 1499]};
end
% if nargin < 3
%     EffectivePixelSize = 0.75;
% end
% if nargin < 4
%     AngleBetweenProjections = 360/1499;
% end
if nargin < 5
    ParentPath = '/mnt/tomoraid-LSDF/rci/MI1057-ESRF-BM05-Sept2011/pc/tomo';
end
%% Set parameters
DataFolderNamePrefix = 'opt_xeno_';
% 0000000001111111111222222222233333333333444444
% 1234567890123456789012345678901234567890123456
% opt_xeno_stage11_tomo_22keV_10x_125mm_0p20sec_
DataFolderStruct = dir([ParentPath '/data/' DataFolderNamePrefix '*']);
%DistanceOffset = 23; %mm
for nn = numel(DataFolderStruct):-1:1
    %NumOfProjs = numel(dir([ParentPath '/int/' DataFolderStruct(nn).name '/int*']));
    %Energy    =  str2double(DataFolderStruct(nn).name(23:24));
    %Distance  = (DistanceOffset + str2double(DataFolderStruct(nn).name(33:35)))/1000;
    EffectivePixelSize = 15/str2double(DataFolderStruct(nn).name(29:30));
    %fprintf('%f %u\n',EffectivePixelSize,NumOfProjs)
    AngleBetweenProjections{nn} = 360/NumOfFirstAndLastProjection{nn}(2);
    %fprintf('%f %u\n',EffectivePixelSize,NumOfFirstAndLastProjection{nn}(2))
    %fprintf('%2u. %-50s: Energy = %2u Distance = %5.3f Pixelsize = %7.2g\n',nn,DataFolderStruct(nn).name,EnergyDistancePixelsize{nn})
end
%% Create .par files.
MakeParFile(StartEndVoxels,NumOfFirstAndLastProjection,EffectivePixelSize,AngleBetweenProjections,ParentPath)