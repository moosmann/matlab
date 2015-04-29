function [binStart count] = ReadHist(filename)

if nargin < 1
    filename = '/mnt/tomoraid-LSDF/tomo/ESRF_MI1079_ID19_July2011_inlineTomo/Xenopus_4cell/vol/Xenopus_4cell_20keV/Histogram of TIElo_alpha5p32.xls';
end

% Open file.
fid = fopen(filename,'r');
% Scan text in the file. The second next line uses a syntax which skips
% the rest of line '%*[^\n]'.
[c] = textscan(fid, '%f %u','Headerlines',1);
fclose(fid);
binStart = c{1};
count = c{2};

% Data set: %s\n'...
%     '%*[^\n]\n%*[^\n]\n%*[^\n]\n%*[^\n]\n%*[^\n]\n'...
%     'Mean detector counts:%f\n'...
%     'Mean variance: %f\n'...
%     'Mean proportionality factor: %f\n'...
%     'Mean real photons: %f\n'...
%     'Mean Poisson noise: %f\n\n
% % Open file to write.
% filename = '/mnt/tomoraid-LSDF/rci/MI1057-ESRF-BM05-Sept2011/pc/tomo/PoissonStatisticsOfFlatFieldsTable.txt';
% fid = fopen(filename,'wt');
% fprintf(fid,'stage E z dx dt <counts> <var> ConvFac <photons> <PoisNoise>\n');
% fprintf(fid,'- keV mm micron ms 1 1 1 1 %%\n');
% DistanceOffset = 13; %mm
% for ff = length(c{1}):-1:1
%     ParStr = char(c{1}(ff));
%     Stage     =  str2double(ParStr(15:16));
%     Energy    =  str2double(ParStr(23:24));%in keV
%     Distance  = (DistanceOffset + str2double(ParStr(33:35)));%in mm
%     Pixelsize = 15/str2double(ParStr(29:30));%in microns
%     ExpoTime  = str2double(ParStr(41:42))*10;%in milliseconds
%     fprintf(1,'%s\n',ParStr);
%     fprintf(fid,'%u %u %u %.2f %u',Stage,Energy,Distance,Pixelsize,ExpoTime);
%     fprintf(fid,' %.1f %.1f %.6f %.1f %.6f\n',c{2}(ff),c{3}(ff),c{4}(ff),c{5}(ff),c{6}(ff));
% end
% fclose(fid);
