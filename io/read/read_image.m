function [im, tif_info] = read_image(filename, filetype, roi, tif_info, shape, dtype)
% Reads 'tif', 'img', 'dar', 'ref', and 'edf' images and all the files
% MATLAB's imread function can read without additional arguments than file
% format.
%
% filename: string. Absolute filename.
% filetype: string. Suffix of image file format. Default: ''
%
% Written by Julian Moosmann. Modified: 2017-06-12
%
% im = read_image(filename, filetype)

%% Default arguments 
if nargin < 2
    filetype = '';
end
if nargin < 3
    roi = [];
end
if nargin < 4
    tif_info = [];
end
if nargin < 5
    shape = [];
end
if nargin < 6
    dtype = '';
end

%% TODO: roi support for edf and tiff files
%% TODO: check flip/rotation/transponation for raw data

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
    case 'tif_ml'
        warning( 'off', 'MATLAB:imagesci:rtifc:missingPhotometricTag');
        switch numel( roi )
            case 0
                im = rot90( imread( filename, 'tif' ), 1);
            case 2
                y0 = max( (3840-roi(2)+1), 1 );
                y1 = min( (3840-roi(1)+1), 3840);
                im = ( rot90( imread( filename, 'tif', 'PixelRegion', {[y0 y1], [1 5120]}), 1) );
            case 4
                y0 = max( (3840-roi(2)+1), 1 );
                y1 = min( (3840-roi(1)+1), 3840);
                x0 = max( (5120-roi(4)+1), 1);
                x1 = min( (5120-roi(1)+1), 5120);
                im = (rot90( imread( filename, 'tif', 'PixelRegion', {[y0 y1], [x0 x1]} ), 1) );
        end
    case 'tif'
        switch numel( roi )
            case 0
                im = flipud( read_tif(filename, tif_info) );
            case 2
                im = flipud( read_tif(filename, tif_info, [(3840-roi(2)+1) (3840-roi(1)+1)] ) );
            otherwise
                error( 'Only vertical ROI supported for fast tif reading.' )
        end
    case {'img', 'dar', 'ref'}        
            im = read_dat_jm( filename, roi );
    case 'raw'
        %% rot90 for FLI but not for KIT
        if shape(1) > shape(2)
            im = flipud( read_raw( filename, shape, dtype, roi ) );
        else            
            %% roi for raw is horizontal roi not verticall!                        
            if isempty( roi )
                im = flipud(rot90(read_raw( filename, shape, dtype, [] ), -1));
            else
                im = read_raw( filename, shape, dtype, [] );
                im = im(roi(1):roi(2),:);
                im = flipud(rot90(im, -1));                
            end
        end        
    otherwise
        im = imread( filename, filetype )';
end
