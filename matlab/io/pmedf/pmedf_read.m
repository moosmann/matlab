function [header, data] = pmedf_read ( f, outputPrecision )
% Reading ESRF data files ('*.edf').
%
% Usage:
%	[header, data] = pmedf_read('hello.edf');
%
% Author: Petr Mikulik
% Version: 11. 5. 2008
% History:
%	    May 2008:
%		Minor clean-up.
%	    June 2006:
%		Use the C++ plugin pmedf_readC if it exists; it is several
%		times faster. (Currently (Octave up to 2.1.72) has slow fread()
%		function for skip=0.)
%	    February 2005:
%		Don't read the data if only the header output argument requested.
%		Better protection against files with a broken header.
%		Added reading .gz and .bz2 files.
%		Reread file if data not fully read (was happening with pipe).
%	    September 2004:
%		if (0) and warning for fscanf('%c') by Octave >=2.1.55.
%	    May 2002:
%		Rewrite into a new routine pmedf_read(); this one does not
%		return a structure of keywords, but the whole header. This may
%		be redefined, later?
%	    2001:
%		pmehf_read updated to read ID19 EDF files.
%	    April 2000:
%		pmehf_read.m for EHF files at ID1.
%       modified by Julian Moosmann, 2013-10-21
%
%[header, data] = pmedf_read ( f, outputPrecision )

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    outputPrecision = 'single';
end
%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if the fast C++ plugin exists, then use it
if exist('pmedf_readC', 'file')
    if nargout == 1
        header = pmedf_readC(f);
    else
        [header, data] = pmedf_readC(f);
    end
    return
end

% open the file
[fid] = fopen(f,'rb');
if fid == -1,
    fprintf('pmedf_read: cannot open file "%s"\n',f);
    return
end

%% read the header
%hsize = 512; % edf header size is a multiple of 512
closing = sprintf('}\n');
% The stupid construct above due to Matlab non-interpretation of escape
% sequences in string; Octave (like C and C++) is OK with
% closing = '}\n';
%[tmp, count] = fread(fid, 512, 'char');
tmp = fread(fid, 512, 'char');
header = sprintf('%c', tmp);
while ~strcmp(header(length(header)-length(closing)+1:length(header)), closing)
    [tmp, count] = fread(fid, 512, 'char');
    header = sprintf('%s%c',header,tmp);
    if count<512 % this is not an edf file
        header = [];
        data = [];
        return;
    end
end

% One output argument requested => return immediately with just the header.
if nargout == 1
    fclose(fid);
    return
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
        dt='double64'; %db=8;
    otherwise
        error(['Unknown data type "', edf.datatype, '" of file "', f, '"']);
end
% Note: dim1 == sizex, dim2 == sizey

% byteorder
switch edf.byteorder
    case 'HighByteFirst'
        arch='ieee-be';
    case 'LowByteFirst'
        arch='ieee-le';
    otherwise
        arch='native';
end

% read binary image
switch lower(outputPrecision)
    case 'double'
        [data, count] = fread(fid, [edf.dim1,edf.dim2], dt, 0, arch);
    case 'single'
        dt = [dt '=>single'];
        [data, count] = fread(fid, [edf.dim1,edf.dim2], dt, 0, arch);
end

if edf.dim1*edf.dim2 ~= count
    % Ooops, a trouble -- data were not completely read. This was sometimes
    % happening when reading many files via pipe.
    % Let's try to reopen and reread the file.
    fprintf('  --> ooops, problem reading last %i B out of %i B, retrying...\n', edf.dim1*edf.dim2-count, edf.dim1*edf.dim2);
    fclose(fid);
    pause(1); % pause 1 second
    %[fid, msg] = fopen(f,'rb');
    fid = fopen(f,'rb');
    if fid == -1,
        fprintf('    --> pmedf_read: cannot reopen file "%s"\n',f);
        return
    end
    % skip the header
    [data, count] = fread(fid, length(header), 'uchar', 0, arch);
    if count~=length(header)
        fprintf('    --> pmedf_read: cannot reread header from file "%s"\n',f);
        return
    end
    % read the data again
    [data, count] = fread(fid, [edf.dim1,edf.dim2], dt, 0, arch);
    if edf.dim1*edf.dim2 ~= count
        fprintf('    --> pmedf_read: cannot reread data from file "%s"\n',f);
        return
    end
end

% data flipping conversion:
% data = data';
% Note: ID19 matlab macros do NOT transpose the image to onze compatibility.

% close file stream
fclose(fid);

% eof pmedf_read.m
