%% Writing ESRF header files: .edf files.
%% Like pmedf_write, it writes the header and the matrix with the data; but
%% in addition, it takes care of the raster axes orientation.
%%
%% Usage:  
%%	new_header = pmedf_writeWithRasterAxes ( filename, header, data, rasterAxes )
%%
%% See 'help pmedf_write' for more details.
%%
%% The addition of axes: by default, pmedf_write() writes the matrix "data" as
%% is, i.e. in matricial system of rasterAxes="XrightYdown", while you would need
%% them in cartesian system of rasterAxes="XrightYup", if the "data" matrix is a
%% cartesian as of data(x,y).
%%
%% Consequently -- examples:
%%	new_header = pmedf_write('bone0007_new.edf', bone.h, bone.a, "XrightYdown");
%%	    ... writes it as matrix
%%	new_header = pmedf_write('bone0007_new.edf', bone.h, bone.a, "XrightYup");
%%	    ... writes it as cartesian "matrix", i.e. flipud'ed.
%%
%% Author: Petr Mikulik
%% Version: 14. 6. 2004
%% History:
%%	June 2004: update -- accept also "matrix" and "cartesian" keywords
%%	June 2003: first version.

function new_header = pmedf_writeWithRasterAxes ( edffile, header, data, rasterAxes )

if nargin ~= 4
  fprintf('Usage: pmedf_writeWithRasterAxes ( filename, header, data, rasterAxes )\n');
  fprintf('where the string rasterAxes can be "XrightYdown" (synonym "matrix") or\n');
  fprintf('"XrightYup" (synonym "cartesian", thus saving a cartesian matrix).');
  return
end

switch rasterAxes
    case 'XrightYdown', ;
    case 'matrix', rasterAxes = 'XrightYdown';
    case 'XrightYup', data = fliplr(data);
    case 'cartesian', data = fliplr(data); rasterAxes = 'XrightYup';
    otherwise error(['Unknown raster axes "', rasterAxes,'"']);
end

header = pmedf_putInHeader(header, 'RasterAxes', rasterAxes);

new_header = pmedf_write(edffile, header, data);

% eof pmedf_writeWithRasterAxes.m
