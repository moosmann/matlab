function out = imgread(filename)
% Read Beckmannsches data format.
%
% Header of 19 bytes following little-endian binary image data.
%
% Arguments
% filename: string. Name of file to read.
% 
% Written by Julian Moosmann. 
% First version: 2016-06-29
% Last modification: 2016-06-29

fid = fopen(filename,'r', 'l');

% Read 19 byte header
header = fscanf(fid, '%c', 19);

% Read dimension from header. TODO: User regexp
dim1 = str2double(header(8:11));
dim2 = str2double(header(13:16));

str = '*uint16';
out = transpose(fread(fid, [dim1, dim2], str));

fclose(fid);