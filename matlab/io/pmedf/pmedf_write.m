%% Writing ESRF header files: .edf files.
%% Writes the header and the matrix with the data. 
%%
%% Usage:  
%%	new_header = pmedf_write ( filename, header, data )
%%
%% Filename should have '.edf', '.edf.gz' or '.edf.bz2' extension.
%%
%% It returns the header of the file written.
%%
%% Data are stored under the short/double/integer format according to
%% the specification in the header.
%%
%% Date key in the header unchanged (yet): Octave is OK, but Matlab does not 
%% provide  "ctime(time());"
%%
%% Example: 
%%	[bone.h, bone.a] = pmedf_read('bone0007.edf');
%%	new_header = pmedf_write('bone0007_new.edf', bone.h, bone.a);
%%
%% Author: Petr Mikulik
%% Version: 5. 6. 2005
%% History:
%%	June 2006: Support for writing signed datatyped.
%%	February 2005: Added writing .gz and .bz2 files.
%%	May 2002: rewrite for string-like header of edf files.
%%	13. 4. 2000: version for ehf files (structure of header fields).

function new_header = pmedf_write ( edffile, header, data )

if nargin ~= 3
  fprintf('Usage: pmedf_write ( filename, header, data )\n');
  return
end

edf.datatype = pmedf_findInHeader(header, 'DataType', 'string');
edf.byteorder = pmedf_findInHeader(header, 'ByteOrder', 'string');
edf.dim1 = pmedf_findInHeader(header, 'Dim_1', 'int');
edf.dim2 = pmedf_findInHeader(header, 'Dim_2', 'int');
edf.size = pmedf_findInHeader(header, 'Size', 'int');

[nr, nc] = size(data);
if edf.dim1~=nr
    header = pmedf_putInHeader(header, 'Dim_1', sprintf('%i',nr));
    edf.dim1 = nr;
end
if edf.dim2~=nc
    header = pmedf_putInHeader(header, 'Dim_2', sprintf('%i',nc));
    edf.dim2 = nc;
end
switch edf.datatype
    case {'UnsignedInteger', 'UnsignedLong'}, dt='uint32'; db=4;
    case 'UnsignedShort', dt='uint16'; db=2;
    case 'UnsignedByte', dt='uint8'; db=1;
    case {'SignedInteger', 'SignedLong', 'Integer'}, dt='int32'; db=4;
    case {'SignedShort', 'Short'}, dt='int16'; db=2;
    case 'SignedByte', dt='int8'; db=1;
    case {'Float', 'FloatValue'}, dt='float32'; db=4;
    case {'Double', 'DoubleValue'}, dt='double64'; db=8;
    otherwise error(['Unknown data type "', edf.datatype, '" of file "', f, '"']);
end
tmp = db*nr*nc; % data size
if edf.size~=tmp
    header = pmedf_putInHeader(header, 'Size', sprintf('%i',tmp));
    edf.size = tmp;
end
switch edf.byteorder
    case 'HighByteFirst', arch='ieee-be';
    case 'LowByteFirst',  arch='ieee-le';
    otherwise arch='native';
end

% open the output file
l = length(edffile);
if l >= 3 & strcmp(edffile(l-2:end),'.gz')	% write .edf.gz
    is_pipe = 1;
    fid = popen(['gzip >', edffile],'w');
elseif l >= 4 & strcmp(edffile(l-3:end),'.bz2') % write .edf.bz2
    is_pipe = 1;
    fid = popen(['bzip2 >', edffile],'w');
else
    is_pipe = 0;
    [fid, msg] = fopen(edffile,'wb');
end
if fid == -1
  fprintf('pmedf_write: cannot write file "%s"\n', edffile);
  return
end

if fid==-1
    fprintf('pmedf_write: cannot write file "%s"\n', edffile);
else
    fprintf('Writing %i x %i x %s to file "%s"\n',edf.dim1,edf.dim2,edf.datatype,edffile);
%    if (ehf.offset ~= 0)
%	fprintf('SKIPPING offset is not yet supported! But it is easy...\n');
%    end
    % Write header:
    fprintf(fid, '%s', header);
    % Write data:
    count = fwrite(fid, data, dt, 0, arch);
%   count = fwrite(fid, data', s, 0, arch);
    if count~=nr*nc
	fprintf('ERROR writing file %s (disk full?)\n', edffile);
    end
    if is_pipe
	pclose(fid);
    else
	fclose(fid);
    end
end

new_header = header;

% eof pmedf_write.m
