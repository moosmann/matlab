%% pmedf_removeInHeader --- remove line of the given key in edf file header
%%
%% Usage:
%%	header = pmedf_removeInHeader( header, key )
%%
%% In the 'header' (a string, usually a multiple of 512 B as the ESRF edf file
%% header), find keyword 'key' and remove the whole associated line. Keep the
%% original size of the header.
%%
%%
%% Tech note: regular expression for a header line is
%%	\nkey [ ]*= [ ]*value[ ]*;\n
%%
%% Example:
%%	header = pmedf_removeInHeader( header, 'row_end' );
%%
%% Author: Petr Mikulik
%% Version: 2. 5. 2002

function header = pmedf_removeInHeader( header, key )

if nargin~=2
    error('syntax: pmedf_removeInHeader(...)');
end

key = ['\n' key ' '];
value = findstr(header, key);
if value==[] return; end

orig_size = length(header);
header_beg = [header(1:value(1)-1) sprintf('\n')];
    % crazy sprintf('\n') instead of '\n' due to Matlab incompatibility
header = header(value(1):length(header));
p = findstr(header, sprintf(';\n'));
header = header(p(1)+2:length(header)-2); % don't copy the trailing '}\n'
header = [header_beg header];

% expand to the original size
header = [header repmat(' ', 1, orig_size-length(header)-2) sprintf('}\n')];

% eof pmedf_removeInHeader.m
