% Open text file. Scan for text pattern. Write scanned text into text file.

% Define filestring.
filestring = '/mnt/tomoraid-LSDF/rci/MI1057-ESRF-BM05-Sept2011/pc/tomo/PoissonStatisticsOfFlatFieldsNoHeading.txt';
% Open file.
fid = fopen(filestring,'r');
% Scan text in the file. The second next line uses a syntax which skips
% the rest of line '%*[^\n]'.
c = textscan(fid, ['Data set: %s\n'...
    '%*[^\n]\n%*[^\n]\n%*[^\n]\n%*[^\n]\n%*[^\n]\n'...
    'Mean detector counts:%f\n'...
    'Mean variance: %f\n'...
    'Mean proportionality factor: %f\n'...
    'Mean real photons: %f\n'...
    'Mean Poisson noise: %f\n\n']);
fclose(fid);
% Open file to write.
filestring = '/mnt/tomoraid-LSDF/rci/MI1057-ESRF-BM05-Sept2011/pc/tomo/PoissonStatisticsOfFlatFieldsTable.txt';
fid = fopen(filestring,'wt');
fprintf(fid,'stage E z dx dt <counts> <var> ConvFac <photons> <PoisNoise>\n');
fprintf(fid,'- keV mm micron ms 1 1 1 1 %%\n');
DistanceOffset = 13; %mm
for ff = length(c{1}):-1:1
    ParStr = char(c{1}(ff));
    Stage     =  str2double(ParStr(15:16));
    Energy    =  str2double(ParStr(23:24));%in keV
    Distance  = (DistanceOffset + str2double(ParStr(33:35)));%in mm
    Pixelsize = 15/str2double(ParStr(29:30));%in microns
    ExpoTime  = str2double(ParStr(41:42))*10;%in milliseconds
    fprintf(1,'%s\n',ParStr);
    fprintf(fid,'%u %u %u %.2f %u',Stage,Energy,Distance,Pixelsize,ExpoTime);
    fprintf(fid,' %.1f %.1f %.6f %.1f %.6f\n',c{2}(ff),c{3}(ff),c{4}(ff),c{5}(ff),c{6}(ff));
end
fclose(fid);
