function a = pmedfread ( fajl )
% Routine  image_matrix = pmedfread(fajl)  is a simplified call to pmedf_read
% which returns two arguments, the image matrix and its header.
%
% Author: Petr Mikulik
% Version: 2002

[~, a] = pmedf_read( fajl );
