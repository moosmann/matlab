function [data] = read_raw( filename, shape, dtype, hor_roi)
% Read arrays of any numerical data type from binary files created with
% read_dat. 
%
% Arguments:
% filename : string, absolute or relative path
% shape : vector, shape of data to read in
% dtype : string. data type
% hor_roi : ROI to read in. Skips the first (hor_roi(1) - 1) lines and reads
% until hor_roi(2)
%
% Written by Julian Moosmann, 2017-11-25

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    filename = '/asap3/petra3/gpfs/p05/2017/data/11003288/raw/syn151_58L_Mg_12_000/syn151_58L_Mg_12_000_ref_0010.raw';
end
if nargin < 2
    shape = [3056 3056];    
end
if nargin < 3
    dtype = 'uint16';
end
if nargin < 4
    hor_roi = [];
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~strcmp( filename(end-2:end), 'raw' )
    filename = sprintf( '%s.raw', filename );
end

% Open file for reading
fid = fopen(filename, 'r');
if (fid < 0)
  error(['Could not open file for reading:'  filename]);
end

% Image shape
si = shape;

% Convert data type to matlab data type
switch dtype
    case 'single' % float32, 'F'        
        byt = 2;
    case 'double' % 'D'        
        byt = 4;
    case 'uint16' %'U'        
        byt = 2;
    case 'L' % 'uint32';        
        byt = 4;
    otherwise
        error(['Data type not supported:' type]);
end
cl = dtype;

if isempty( hor_roi )
    [data, cnt] = fread(fid, si, [cl '=>' cl]);
    if prod(si) ~= cnt
        error('Number of elements of data read mismatches provided image shape.')
    end
else
    si = [si(1), (hor_roi(2) - hor_roi(1) + 1)];
    fseek( fid, ( hor_roi(1) - 1 ) * si(1) * byt, 0);
    [data, cnt] = fread(fid, si, [cl '=>' cl]);
    if prod(si) ~= cnt
        error('Number of elements of data read mismatches provided image shape.')        
    end
end

% Close the file
fclose(fid);

end
