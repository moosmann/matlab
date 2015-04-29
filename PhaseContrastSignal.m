function out = PhaseContrastSignal(NumProjAv,doShow)
% Calculate mean value and standard deviation of a central region of the
% xenopus data and print it to a file. The mean values and standard
% deviations are averages of several projections.

%% Default arguments.
if nargin < 1
    NumProjAv = 20;
end
if nargin < 2
    doShow = 0;
end

%% Loop over data sets found under ParentPath.
ParentPath = '/mnt/tomoraid-LSDF/rci/MI1057-ESRF-BM05-Sept2011/pc/tomo/';
DataFolderNamePrefix = 'opt_xeno_';
DataFolderStruct = dir([ParentPath 'int/' DataFolderNamePrefix '*']);
NumDataSets = numel(DataFolderStruct);
out(NumDataSets).mean = 0;
out(NumDataSets).std = 0;
DistanceOffset = 13; %mm
% 1/CropFac*Dim will be cropped center symmetrically.
CropFac = 3;
fprintf('\n')
% Define text file name and open file to write.
txtFileName = [ParentPath 'PhaseContrastSignal2.txt'];
fprintf('Text file: %s\n',txtFileName)
fid = fopen(txtFileName,'wt');
fprintf(fid,'stage E z dx dt mean std std*mean std2mean min max max-min\n');
fprintf(fid,'- keV mm micron ms 1 1 1 1 1 1 1\n');
for nn = 1:NumDataSets
    % 0000000001111111111222222222233333333333444444
    % 1234567890123456789012345678901234567890123456
    % opt_xeno_stage11_tomo_22keV_10x_125mm_0p20sec_
    DataFolderName = DataFolderStruct(nn).name;
    DataPath       = [ParentPath 'int/' DataFolderName '/'];
    %% Read parameters from directory names.
    out(nn).Stage     =  str2double(DataFolderName(15:16));
    out(nn).Energy    =  str2double(DataFolderName(23:24));%in keV
    out(nn).Distance  = (DistanceOffset + str2double(DataFolderName(33:35)));%in mm
    out(nn).Pixelsize = 15/str2double(DataFolderName(29:30));%in microns
    out(nn).ExpoTime  = str2double(DataFolderName(41:42))*10;%in milliseconds
    %fprintf('%2u. %-50s: ',nn,DataFolderName);
    fprintf('Stage:%2u Energy:%2u Pixelsize:%4.02g Distance:%3u ExposureTime:%3u',...
        out(nn).Stage,out(nn).Energy,out(nn).Pixelsize,out(nn).Distance,out(nn).ExpoTime)
   
    %% Calculate mean values and standard deviations.
    fs       = dir([DataPath 'int*.edf']);
    NumFiles = length(fs)-6;
    for ff = NumProjAv:-1:0
        % Read image and get dimension.
        im = pmedfread([DataPath fs(floor(NumFiles/NumProjAv*ff)+1).name])';
        [dimx dimy] = size(im);
        % Crop central region.
        xx = ceil(dimx/2*(CropFac-1)/CropFac):floor(dimx/2*(CropFac+1)/CropFac);
        yy = ceil(dimy/2*(CropFac-1)/CropFac):floor(dimy/2*(CropFac+1)/CropFac);
        im = im(xx,yy);
        imMean(ff+1) = mean(im(:));
        imStd(ff+1)  = std(im(:));
        imMax(ff+1)  = max(im(:));
        imMin(ff+1)  = min(im(:));
    end
    % Average of mean values and standard deviations.
    out(nn).mean     = mean(imMean(:));
    out(nn).std      = mean(imStd(:));
    out(nn).stdmean  = out(nn).std*out(nn).mean;
    out(nn).std2mean = out(nn).std/out(nn).mean;
    out(nn).max      = mean(imMax(:));
    out(nn).min      = mean(imMin(:));
    fprintf(fid,'%u %u %u %.2f %u %g %g %g %g %g %g %g\n',...
        out(nn).Stage,out(nn).Energy,out(nn).Distance,out(nn).Pixelsize,out(nn).ExpoTime,out(nn).mean,out(nn).std,out(nn).stdmean,out(nn).std2mean,out(nn).min,out(nn).max,out(nn).max-out(nn).min);
    % Print more info.
    if doShow
        ishow(im)
        fprintf('MEANs: %s\n',mat2str(imMean,3))
        fprintf('STDs: %s\n',mat2str(imStd,3))
    end
    fprintf(' MEAN:%5.3g',out(nn).mean)
    fprintf(' STD:%6.3g',out(nn).std) 
    fprintf(' STD*MEAN:%7.3g',out(nn).stdmean)
    fprintf(' STD/MEAN:%6.3g',out(nn).std2mean)
    fprintf(' MIN:%5.3g',out(nn).min) 
    fprintf(' MAX:%5.3g',out(nn).max)
    fprintf(' MAX-MIN:%5.3g\n',out(nn).max-out(nn).min)
end
fclose(fid);
