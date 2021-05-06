function [data] = read_dat_jm( filename, vert_roi)
% Read arrays of any numerical data type from binary files created with
% read_dat. 
%
% Arguments:
% filename : string, absolute or relative path
% vert_roi : ROI to read in. Skips the first (vert_roi(1) - 1) lines and reads
% until vert_roi(2)
%
% Written by Tilman Donath, 7.1.2007. Modified by JH, 11.12.09. Tilman
% Donath debugged the file selection, 22.1.2009.
% Rewritten by Julian Moosmann, 2017-05-20

%% TODO: Test vert_roi and adopt for data other than 2D.

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    [filename, pathname, ~] = uigetfile({'*.dat'; '*.*'});
    if filename(1) == 0
      fprintf('No valid file selected. Returning.')
      return
    end
  filename = [pathname filename];
  fprintf(['Selected file: \n ' filename '\n']);
end
if nargin < 2
    vert_roi = [];
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Open file for reading
%fid = fopen(filename, 'r');
permission = 'r';
machinefmt = 'n';
encodingIn = 'UTF-8';
fid = fopen(filename,permission,machinefmt,encodingIn);
if (fid < 0)
  error(['Could not open file for reading:'  filename]);
end

% Read header and interpret
header = fgets(fid);

pos = strfind(header, '_');
if (numel(pos) < 3) 
    error('Wrong header information.');
end

% Image dimension and type
ndim = str2double(header((pos(1) + 1):(pos(2) - 1)));
type = header((pos(2) + 1):(pos(3) - 1));

if (numel(pos) < (3 + ndim)) 
    error('Wrong header information.');
end

% Image shape
si = zeros(1,2);
for i_dim = 0:(ndim - 1)    
    si(i_dim+1) = str2double(header((pos(3 + i_dim) + 1):(pos(3 + i_dim + 1) - 1)));
end

switch ndim
    case 0
         si =[1 1];
    case 1
         si =[1 si];
end

% Convert data type to matlab data type
switch type
    case 'F'
        cl = 'single'; %(float32)
        byt = 2;
    case 'D'
        cl = 'double';
        byt = 4;
    case 'U'
        cl = 'uint16';
        byt = 2;
    case 'L'
        cl = 'uint32';
        byt = 4;
    case 'I'
        byt = 'int16';
    otherwise
        error(['Data type not supported:' type]);
end

if isempty( vert_roi )
    [data, cnt] = fread(fid, si, [cl '=>' cl]);
    if prod(si) ~= cnt
        error('Number of elements in read data mismatches header information.')
    end
else
    si = [si(1), (vert_roi(2) - vert_roi(1) + 1)];
    fseek( fid, ( vert_roi(1) - 1 ) * si(1) * byt, 0);
    [data, cnt] = fread(fid, si, [cl '=>' cl]);
    if prod(si) ~= cnt
        error('Number of elements in read data mismatches header information.')  
    end
end

% Close the file
fclose(fid);

end
