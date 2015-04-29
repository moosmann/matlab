%% Make a diff series of a series of edf files. It saves the first file as it
%% is, and then it saves only differences with respect to the previous image.
%% It can compress output files on fly by gzip or bzip2.
%% 
%% Results: you should test whether it saves space or not on real edf data
%% (radiography, diffraction images, etc.); then write the reverse summing
%% procedure. Actually, whenever I tried to use this procedure, it did not
%% make the resulting compressed file smaller. That's probably because there
%% is always some noise and thus there are no large ranges of zeros in the
%% difference image; and the bzip2 compression algorithm is pretty good!
%%
%% Syntax:
%%	pmedf_diff_series( outprefix, outsuffix, edfformat, numbers )
%% where
%%	outprefix - can be a directory name (with slash)
%%	outsuffix - '' for saving full edf files, or '.gz' or '.bz2' for
%%		    compressing them on fly
%%	edfformat - format string for input file names
%%	numbers	  - range of file name numbers
%%
%% Examples:
%%	pmedf_diff_series( 'new/', '',     'tomo%03i.edf', [0:800] );
%%	pmedf_diff_series('t/',    '.bz2', '/data/scan01/x%03i.edf', [0:4]);
%%	pmedf_diff_series('t/new', '.gz',  'x%03i.edf', [0:2:1200]);
%%
%% Note: currently hardcoded header size of 1024 B.
%%
%% Author: Petr Mikulik
%% Version: May 2002
%% History:
%%    1. 5. 2002: Updated for reading files from not current directory.
%%   27. 4. 2002: First version.

function pmedf_diff_series( outprefix, outsuffix, edfformat, numbers )

if nargin ~= 4
  fprintf ('Usage: pmedf_diff_series( outprefix, outsuffix, edfformat, numbers )\n');
  return
end

% *** Setup - hardcoded limits ***
headerlen = 1024;

% Read the first edf file
name = sprintf(edfformat, numbers(1));
[fid,msg] = fopen(name,'rb');
if fid == -1,
    fprintf('pmedf_diff_series: cannot open file "%s"\n',name);
    return
end

% read the header of the 1st file
[header, count] = fread(fid, headerlen, 'char', 0);
header = sprintf('%c',header);
fclose(fid);

if     index(header,'UnsignedInteger') dtype = 'uint32';
elseif index(header,'UnsignedLong') dtype = 'uint32';
elseif index(header,'UnsignedShort') dtype = 'uint16';
elseif index(header,'UnsignedByte') dtype = 'uint8';
elseif index(header,'Float') | strcmp(header,'FloatValue') dtype = 'float';
else dtype = 'x';
     fprintf('Unknown binary data type in the header!\n');
end

if     index(header,'HighByteFirst') arch = 'ieee-be';
elseif index(header,'LowByteFirst')  arch = 'ieee-le';
else   arch = 'native';
end

% now proceed over all files
for k=1:length(numbers)
    name = sprintf(edfformat, numbers(k));
    nodirname = rindex(name, '/');
    nodirname = name(nodirname+1:length(name));
    fprintf('Input file %s\t',name);
    [fid,msg] = fopen(name, 'rb');
    if fid == -1,
	fprintf('pmedf_diff: cannot open file "%s"\n',name);
    end
    % read header of the file
    [header, count] = fread(fid, headerlen, 'char', 0);
%   header = sprintf('%c',header);
    % remember previous data
    if k>1
	prev_data = data;
    end
    % read the image
    [data, count] = fread(fid, Inf, dtype, 0, arch);
    fclose(fid);
%    data = data';  % for presentation purposes; also may reshape...

    % diff data
    if k>1
	diff_data = data - prev_data;
    else
	diff_data = data;
    end

    % write the output file; firstly, determine the name and open it
    outname = sprintf('%s%s%s', outprefix, nodirname, outsuffix);
    if strcmp(outsuffix,'')
	[fid,msg] = fopen( outname, 'wb' );
	if fid == -1
	    fprintf('Cannot write file "%s (disk full? protected?)"\n', outname);
	    return
	end
    elseif strcmp(outsuffix,'.bz2')
	outname = sprintf('bzip2 -c >%s', outname);
	fid = popen( outname, 'w' ); % supposing binary mode
    elseif strcmp(outsuffix,'.gz')
	outname = sprintf('gzip -c >%s', outname);
	fid = popen( outname, 'w' ); % supposing binary mode
    else
	error( sprintf('Unknown suffix %\n', outsuffix) );
    end
    fprintf('=> output to %s\n', outname);
    if fid == -1
	fprintf('Cannot write file "%s (disk full? broken pipe? protected?)"\n', outname);
	return
    end
    % write data	
    count = fwrite(fid, header, 'char', 0);
    count = fwrite(fid, diff_data, dtype, 0, arch);
    fclose(fid);
end


% endfunction

%eof pmedf_diff_series.m
