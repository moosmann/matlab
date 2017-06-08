function im = read_image(filename, filetype, roi)
% Reads 'tif', 'img', 'dar', 'ref', and 'edf' images and all the files
% MATLAB's imread function can read without additional arguments than file
% format.
%
% filename: string. Absolute filename.
% filetype: string. Suffix of image file format. Default: ''
%
% Written by Julian Moosmann. Modified: 2016-11-07
%
% im = read_image(filename, filetype)

%% Default arguments 
if nargin < 2
    filetype = '';
end
if nargin < 3
    roi = [];
end

%% TODO: roi support for edf and tiff files

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty( filetype )
    [~, ~, filetype] = fileparts(filename);
end
if strcmp( filetype(1), '.')
    filetype(1) = [];
end

switch lower( filetype )
    case 'edf'
        im = pmedfread( filename );
    case {'tif', 'tiff'}
        warning( 'off', 'MATLAB:imagesci:rtifc:missingPhotometricTag');
        switch numel( roi )
            case 0
                im = rot90( imread( filename, 'tif' ), -1);
            case 2
                y0 = max( (3840-roi(2)), 0 );
                y1 = min( (3840-roi(1)), 3840);
                im = rot90( imread( filename, 'tif', 'PixelRegion', {[y0 y1], [1 5120]}), -1);
            case 4
                y0 = max( (3840-roi(2)), 0 );
                y1 = min( (3840-roi(1)), 3840);
                x0 = max( (5120-roi(4)), 0);
                x1 = min( (5120-roi(1)), 5120);
                im = rot90( imread( filename, 'tif', 'PixelRegion', {[y0 (3840-roi(1))], [x0 (5120-roi(3))]} ), -1);
        end
    case {'img', 'dar', 'ref'}        
            im = read_dat_jm( filename, roi );        
    otherwise
        im = imread( filename, filetype )';
end
