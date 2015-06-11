function [header, a] = pmedf_average4 ( header, a )
% pmedf_average4 - average over 2x2 pixels of the matrix, and remove (later
% versions may update) the corresponding lines in the edf header.
%
% Usage:
%	[header, data] = pmedf_average4 ( header, data );
%
% Author: Petr Mikulik
% Version: 1. 5. 2002
% History:
%	 1. 5. 2002: pmedf_average4.m - rewrite for edf files.
%	13. 4. 2000: average4ehf.m - averaging ehf files.

if nargin ~= 2
  error('usage: pmedf_average4(...)');
  return
end

nrc = floor(size(a)/2);
a = a([1:2:2*nrc(1)],:) + a([2:2:2*nrc(1)],:); % sum odd and even rows
a = a(:,[1:2:2*nrc(2)]) + a(:,[2:2:2*nrc(2)]); % sum odd and even columns
a = a*0.25;

header = pmedf_removeInHeader(header,'row_beg');
header = pmedf_removeInHeader(header,'row_end');
header = pmedf_removeInHeader(header,'col_beg');
header = pmedf_removeInHeader(header,'col_end');
rbin = pmedf_findInHeader(header, 'row_binning', 'int');
if rbin==[] rbin = 1; end
header = pmedf_putInHeader(header,'row_binning', sprintf('%i', 2*rbin), 16);
cbin = pmedf_findInHeader(header, 'col_binning', 'int');
if cbin==[] cbin = 1; end
header = pmedf_putInHeader(header,'col_binning', sprintf('%i', 2*cbin), 16);
% note: optics_used is unchanged
% pixel sizes: unchanged now
%   ehf.psize_x = 2*ehf.psize_y;	% pixel size is double
%   ehf.psize_y = 2*ehf.psize_y;
% scale_x, scale_y: unchanged now
% center_x, center_y: unchanged now

% eof pmedf_average4.m
