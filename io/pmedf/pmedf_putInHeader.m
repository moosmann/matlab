%% pmedf_putInHeader --- puts the given key and its value in edf file header:
%% a new entry is added at the end, an existing key is updated.
%%
%% Usage:
%%	header = pmedf_putInHeader( header, key, new_value [, pos] )
%%
%% In the 'header' (a string, usually a multiple of 512 B as the ESRF edf file
%% header), find keyword 'key'; if it exists, then replace its value by
%% 'new_value', otherwise add a new line. Keep the original size of the
%% header if possible, otherwise extend it by multiples of 512 B.
%% The "=" separating the key and its value is copied from the existing
%% position if the key exists, otherwise positioned at the given 'pos'
%% or separated from the key by a single space.
%%
%% Tech note: regular expression for a header line is
%%	\nkey [ ]*= [ ]*value[ ]*;\n
%%
%% Examples:
%%	header = pmedf_putInHeader( header, 'energy', 12.3 );
%%	header = pmedf_putInHeader( header, 'note', 'after alignment' );
%%	header = pmedf_putInHeader( header, 'note', 'after alignment', 16 );
%%
%% Author: Petr Mikulik
%% Version: 23. 9. 2004
%% History: September 2004: Use isempty() where necessary.
%%	    May 2002: First version.

function header = pmedf_putInHeader( header, key, new_value, pos )

if nargin<3 | nargin>4
    error('syntax: pmedf_putInHeader(...)');
end

key = [sprintf('\n') key ' '];
    % bloody Matlab \n incompatibility
if isnumeric(new_value) new_value = sprintf('%g', new_value); end
value = findstr(header, key);
orig_size = length(header);

if ~isempty(value) % the key already exists
    header_beg = [header(1:value(1)-1)];
    header = header(value(1):length(header));
    p = findstr(header, sprintf(';\n'));
    line = header(1:p(1)-1);
    header = header(p(1)+2:length(header)-2); % don't copy the trailing '}\n'
    p = findstr(line, '= ');
    line = line(1:p(1)+1);
    header = [header_beg line new_value sprintf(' ;\n') header];
    % strip the trailing spaces after the last line
    p = rindex(header, sprintf(';\n'));
    header = header(1:p+1);
else % the key does not exist yet
    p = rindex(header, sprintf(';\n'));
    header = header(1:p);
    if nargin==3 
	pos = 0;
    end
    pos = pos - length(key) - 1; % how many spaces to add
    if pos > 1
	spaces = [repmat(' ', 1, pos+1)];
	header = [header key spaces '= ' new_value sprintf(' ;\n')];
    else
	header = [header key '= ' new_value sprintf(' ;\n')];
    end
end

% now extend the header to its original size or extend it
l = length(header);
if l+2 > orig_size % extend the header
    orig_size = 512*ceil(l/512);
end
header = [header repmat(' ', 1, orig_size-l-2) sprintf('}\n')];

% eof pmedf_putInHeader.m
