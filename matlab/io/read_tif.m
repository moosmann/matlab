function [im, tiff_info] = read_tif( filename, tiff_info, vert_roi )
% Read tiff image using low-level file read. Twice as fast MATLAB's imread.
%
% ARG:
% filename : str
% tiff_info : tiff info struct provided by imfinfo
% vert_roi : 2-vector of integers defining the vertical ROI
%
% Written by J. Moosmann

%% Defaults %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    filename = '/asap3/petra3/gpfs/p05/2017/data/11002839/raw/ehh_2017_019_f/proj_0000.tif';
end
if nargin < 2
    tiff_info = imfinfo( filename, 'tif' ) ;
end
if nargin < 3
    vert_roi = [];
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty( tiff_info )
    tiff_info = imfinfo( filename, 'tif' ) ;
end

% Open file for reading
fid = fopen( filename, 'r');
if (fid < 0)
  error(['Could not open file for reading:'  filename]);
end

switch tiff_info.BitsPerSample    
    case 16
        cl = 'uint16';
        byt = 2;
    case 32
        cl = 'uint32';
        byt = 4;
    otherwise
        error(['Data type not supported:' type]);
end

im_shape = [tiff_info.Width, tiff_info.Height];
offset = tiff_info.StripOffsets(1);

% Read data
if isempty( vert_roi )
    fseek( fid, offset, 0);
    [im, cnt] = fread(fid, im_shape, [cl '=>' cl]);
    if prod(im_shape) ~= cnt
        error('Number of elements in read data mismatches header information.')
    end
else
    fseek( fid, offset, 0);
    im_shape = [im_shape(1), (vert_roi(2) - vert_roi(1) + 1)];
    fseek( fid, ( vert_roi(1) - 1 ) * im_shape(1) * byt, 0);
    [im, cnt] = fread(fid, im_shape, [cl '=>' cl]);
    if prod(im_shape) ~= cnt
        error('Number of elements in read data mismatches header information.')  
    end
end

% Close the file
fclose(fid);

end
