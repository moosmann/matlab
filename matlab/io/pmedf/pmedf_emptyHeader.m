%% Returns an empty header (1x1x1 B).
%%
%% Syntax:
%%	pmedf_emptyHeader ( [headersize] );
%% where
%%	headersize: an optional parameter (default 1024) extends the header
%%	size to this value; note this it has to be a multiple of 512.
%%	Note: special case were headersize==0 can be used for writing EHF
%%	files (separate header and binary data file); then, no padding is
%%	performed. Use with caution.
%%
%%
%% Author: Petr Mikulik
%% Version: May 2002

function header = pmedf_emptyHeader ( headersize )

if nargin > 1
    error('usage: pmedf_emptyHeader( [headersize] )');
end

if nargin==0
    headersize = 1024;
end

header = sprintf('{\nHeaderID = EH:000001:000000:000000 ;\nImage = 1 ;\nByteOrder = LowByteFirst ;\nDataType = UnsignedByte ;\nDim_1 = 1 ;\nDim_2 = 1 ;\nSize = 1 ;\n');
% sprintf('\n') due to Matlab incompatibility

if headersize==0
    header = [header sprintf('}\n')];
else
    if mod(headersize, 512)~=0
	headersize = 512*ceil(headersize/512);
    end
    if length(header)+2 > headersize
	headersize = 512*ceil((length(header)+2)/512);
    end
    header = [header repmat(' ', 1, headersize-length(header)-2) sprintf('}\n')];
end

% eof pmedf_emptyHeader.m
