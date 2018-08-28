function [data] = read_binary_data( filename, shape, dtype, machinefmt)
% Read arrays of any numerical data type from binary files created with
% read_dat.
%
% Arguments:
% filename : string, absolute or relative path
% shape : vector, shape of data to read in
% dtype : string. data type
%
% Written by Julian Moosmann, 2017-11-25

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    filename = '/asap3/petra3/gpfs/p05/2017/data/11003288/processed/syn166_104R_Mg10Gd_4w/optical_flow/syn166_results_to_send/out_1/vec_x_1_1240x1120x20_32bit.vol';
end
if nargin < 2
    shape = [1240 1120 20];
end
if nargin < 3
    dtype = 'single';
end
if nargin < 4
    machinefmt = 'l';
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if ~strcmp( filename(end-2:end), 'raw' ) && strcmp( filename(end-2:end), 'vol' )
%     filename = sprintf( '%s.raw', filename );
% end

% Open file for reading
fid = fopen(filename, 'r');
if (fid < 0)
    error(['Could not open file for reading:'  filename]);
end

% Number of bytes required in the case of ROI reading implementation
%byte = bytes_of_dtype( dtype );

if numel( shape ) == 2
    % Read
    [data, cnt] = fread(fid, shape, [dtype '=>' dtype], 0, machinefmt);
else
    % Read
    [data, cnt] = fread(fid, inf, [dtype '=>' dtype], 0, machinefmt);
    
    % Reshape
    data = reshape( data, shape );
end

% Close the file
fclose(fid);

if prod(shape) ~= cnt
    error('Number of elements of data read mismatches provided image shape.')
end

end
