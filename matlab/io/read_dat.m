% Name: read_dat
%
% Description:
% Read arrays of any numerical data type from binary files created with
% read_dat. 
%
% Parameters:
%
% History: Written by Tilman Donath, 7.1.2007 
%
%  22.1.2009, Tilman Donath: Debugged the file selection
% Modified by Julian Moosmann, 2017-05-20


function [data] = read_dat(filename)

if nargin < 1
    [filename, pathname, ~] = uigetfile({'*.dat'; '*.*'});
    if filename(1) == 0
      fprintf('No valid file selected. Returning.')
      return
    end
  filename = [pathname filename];
  fprintf(['Selected file: \n ' filename '\n']);
end


% Open file for reading
%fid = fopen(filename, 'r');
permission = 'r';
machinefmt = 'n';
encodingIn = 'UTF-8';
fid = fopen(filename,permission,machinefmt,encodingIn);
if (fid < 0)
  error(['Could not open file for reading:'  filename]);
end

% Read header and interprete
header = fgets(fid);

pos = strfind(header, '_');
if (numel(pos) < 3) 
    error('Wrong header information.');
end

% comp = header(1:(pos(1) - 1));
ndim = str2double(header((pos(1) + 1):(pos(2) - 1)));
type = header((pos(2) + 1):(pos(3) - 1));

if (numel(pos) < (3 + ndim)) 
    error('Wrong header information.');
end

% shape=[];
% for i_dim = 0:(ndim - 1)
%     thisdim = str2double(header((pos(3 + i_dim) + 1):(pos(3 + i_dim + 1) - 1)));
%     shape = [shape, thisdim];
% end
shape = [];
for i_dim = (ndim - 1):-1:0
    shape(i_dim+1) = str2double(header((pos(3 + i_dim) + 1):(pos(3 + i_dim + 1) - 1)));
end

% Header information by user included
% if (numel(pos) > (3 + ndim))
%     header_info = header((pos(3 + ndim) + 1):(numel(header) - 1));
% end


% Convert data type to matlab data type
switch type
    case 'F'
        dtype = 'single'; %(float32)
    case 'D'
        dtype = 'double';
    case 'U'                %corrected for HARWI-data type
        dtype = 'uint16';      %JH, 11.12.09
    case 'L'
        dtype = 'uint32';
    case 'I'
        dtype = 'int16';
    otherwise
        error(['Data type not supported:' type]);
end

% Read the data in binary format
%fwrite(fid, data, dtype);
%[data, count] = fread(fid, inf, dtype);
%STATUS = fseek(fid, 1, 'cof') 
[data, cnt] = fread(fid, inf, [dtype '=>' dtype]);

% Close the file
fclose(fid);

% Reshape the data, for Matlab force two-dimensional array
if (ndim == 0) 
    shape =[1 1];
end
if (ndim == 1)
    shape =[1 shape];
end
%if (prod(shape) ~= numel(data))
if prod(shape) ~= cnt
    error('Number of elements in read data mismatches header information.')
end

data = reshape( data, shape);

end
