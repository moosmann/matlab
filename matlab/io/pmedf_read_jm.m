function [data, header] = pmedf_read_jm (f,outputDataType )
% Reading ESRF header files: .edf files.
%
% Usage:
%	[data header] = pmedf_read('hello.edf');

%		Use the C++ plugin pmedf_readC if it exists; it is several
%		times faster. (Currently (Octave up to 2.1.72) has slow fread()
%		function for skip=0.)
% Note: dim1 == sizex, dim2 == sizey
%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    fprintf('Usage:\n');
    fprintf('  [header {, image}] = pmedf_read.m(edf_filename)\n');
    return
end
if nargin < 2
    outputDataType = 'double';
end
%% Body %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% open the file
[fid] = fopen(f,'rb');
if fid == -1,
    fprintf('pmedf_read: cannot open file "%s"\n',f);
    return
end

%% Read the header info
%hsize = 512; % edf header size is a multiple of 512
closing = sprintf('}\n');
% The stupid construct above due to Matlab non-interpretation of escape
% sequences in string; Octave (like C and C++) is OK with
% closing = '}\n';
tmp = fread(fid, 512, 'char');
header = sprintf('%c', tmp);
while ~strcmp(header(length(header)-length(closing)+1:length(header)), closing)
    [tmp, count] = fread(fid, 512, 'char');
    header = [header, sprintf('%c', tmp)];
    if count<512 % this is not an edf file
        header = [];
        data = [];
        return;
    end
end

edf.datatype = pmedf_findInHeader(header, 'DataType', 'string');
edf.byteorder = pmedf_findInHeader(header, 'ByteOrder', 'string');
edf.dim1 = pmedf_findInHeader(header, 'Dim_1', 'int');
edf.dim2 = pmedf_findInHeader(header, 'Dim_2', 'int');
edf.size = pmedf_findInHeader(header, 'Size', 'int');

% Read the binary file
switch edf.datatype
    case {'UnsignedInteger', 'UnsignedLong'}
        dt='uint32'; %db=4;
    case 'UnsignedShort'
        dt='uint16'; %db=2;
    case 'UnsignedByte'
        dt='uint8'; %db=1;
    case {'SignedInteger', 'SignedLong', 'Integer'}
        dt='int32'; %db=4;
    case {'SignedShort', 'Short'}
        dt='int16'; %db=2;
    case 'SignedByte'
        dt='int8'; %db=1;
    case {'Float', 'FloatValue'}
        dt='float32'; %db=4;
    case {'Double', 'DoubleValue'}
        dt='float64'; %db=8;
    otherwise
        error(['Unknown data type "', edf.datatype, '" of file "', f, '"']);
end

%fprintf('Reading %i x %i x %s from file \"%s\"\n',edf.dim1,edf.dim2,dt,f);

%% Byteorder
switch edf.byteorder
    case 'HighByteFirst'
        arch='ieee-be';
    case 'LowByteFirst'
        arch='ieee-le';
    otherwise
        arch='native';
end

data = fread(fid, [edf.dim1,edf.dim2],[dt '=>' outputDataType], 0, arch);

fclose(fid);
% eof pmedf_read_jm.m
