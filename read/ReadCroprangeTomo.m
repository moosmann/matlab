function [out c] = ReadCroprangeTomo(PathToCroprangeFile)
% Check if a croprange file exist which defines the first and last pixels of
% the area that should be cropped. If it not exist the output is 0. If it
% exist the output is 2x2-matrix as follows out = [[HorFirst
% HorLast];[VerFirst VerLast]]; c is a cell and output for debugging.

%% Default arguments.
if nargin < 1
    PathToCroprangeFile = '/mnt/tomoraid-LSDF/tomo/APS_2BM_LifeCellImaging_GUP28266/lifecellimaging/vol/wildtype_30keV_10min_deadtime_25tomo_stage11p0_upwards_620mm_015ms';
end
if nargin < 2
    CroprangeFileName = 'croprange.txt';
end
%% Body.
if PathToCroprangeFile(end) ~= '/'
    PathToCroprangeFile = [PathToCroprangeFile '/'];
end
FilePath = [PathToCroprangeFile CroprangeFileName];
% Check if a croprange file is existing.
if exist(FilePath,'file')
    % Open file.
    fid = fopen(FilePath,'r');
    % Example of croprange file:
    % dir\pixel   first   last    offset
    % horizontal  13      852     0
    % vertical    163     852     0
    c = textscan(fid,'%*s %*s %*s %*s\n%*s %u %u %u\n%*s %u %u %u\n');
    % Close file.
    fclose(fid);
    % HorFirst  = c{1};
    % HorLast   = c{2};
    % HorOffset = c{3};
    % VerFirst  = c{4};
    % VerLast   = c{5};
    % VerOffset = c{6};
    if isempty(c{6})
        c{6} = 0;
    end
    % Define output matrix.
    out = double([[c{1} c{2} c{3}];[c{4} c{5} c{6}]]);
else
    out = 0;
    c = 0;
end
