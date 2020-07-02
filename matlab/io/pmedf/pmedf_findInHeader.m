%% pmedf_findInHeader --- find value of a given key in edf file header
%%
%% Usage:
%%	value = pmedf_findInHeader( header, key [, typ] )
%%
%% In the 'header' (a string, usually a multiple of 512 B as the ESRF edf file
%% header), find keyword 'key' and return its 'value' converted to type 'typ',
%% which can be 'string' (default), 'int' or 'float'.
%%
%% Returns [] if the keyword has not been found.
%%
%% Tech note: regular expression for a header line is
%%	\nkey [ ]*= [ ]*value[ ]*;\n
%%
%% Examples:
%%	pmedf_findInHeader( header, 'Dim_1', 'int' );
%%	pmedf_findInHeader( header, 'Title', 'string' );
%%	pmedf_findInHeader( header, 'Title' );
%%	pmedf_findInHeader( header, 'PSize_1', 'float' );
%%
%% Author: Petr Mikulik
%% Version: 23. 9. 2004
%% History: September 2004: Use isempty() where necessary.
%%	    May 2002: First version.

function value = pmedf_findInHeader( header, key, typ )

if nargin<2 | nargin>3
    error('syntax: pmedf_findInHeader(...)');
end

key = sprintf(['\n' key ' ']);
% again, the above sprintf() due to bloody Matlab \n-interpretation
value = findstr(header, key);
if isempty(value) return; end

header = header(value(1):length(header));
p = findstr(header, sprintf(';\n'));
header = header(1:p(1)-1);
p = findstr(header, '= ');
header = header(p(1)+2:length(header));
% trim trailing spaces
value = deblank(header);
% trim leading spaces
p = findstr(value, ' ');
if ~isempty(p)
    % Octave OK, Matlab fails: p = find(p - [1:length(p)] > 0)(1);
    p = find(p - [1:length(p)] > 0);
    value = value(p(1):length(value));
end

% Return value according to 'typ':
if nargin==2 | strcmp(typ,'string') return; end
if strcmp(typ,'int') [value, tmp] = sscanf(value, '%i', 1); return; end
if strcmp(typ,'float') [value, tmp] = sscanf(value, '%g', 1); return; end
error('unknown "typ"');

% eof pmedf_findInHeader.m
