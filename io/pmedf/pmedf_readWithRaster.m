%% Reading ESRF header files: .ehf / .edf files, and take care of the data
%% matrix orientation.
%%
%% It reads an edf file as pmedf_read() does, but then looks to the header
%% field 'RasterAxes', and if it finds 'XrightYup' or 'cartesian', then does
%% the necessary fliplr(data).
%%
%% Usage:
%%	[header, data] = pmedf_readWithRaster('hello.edf');
%%
%% Author: Petr Mikulik
%% Version: 14. 6. 2004

function [header, data] = pmedf_readWithRaster( f )

if nargin ~= 1
  fprintf ('Usage:\n');
  fprintf ('  [header, image] = pmedf_readWithRaster.m(edf_filename)\n');
  return
end

[header, data] = pmedf_read(f);

r = pmedf_findInHeader( header, 'RasterAxes', 'string' );

if (strcmp(r,'XrightYup') || strcmp(r,'cartesian'))
    data = fliplr(data);
end

%eof pmedf_readWithRaster.m
