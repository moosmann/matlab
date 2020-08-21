function [im, tif_info] = read_image(filename, filetype, roi, tif_info, shape, dtype, trafo, tifftrafo)
% Read 'tif', 'img', 'dar', 'ref', and 'edf' images and all formats that
% are MATLAB's imread function can read without additional arguments than file
% format.
%
% filename: string. Absolute filename.
% filetype: string. Suffix of image file format. Default: ''
% trafo : string, default ''. string which will be evaluated on the read
% image, e.g. 'rot90(im)'. Note that 'im' is the mandatory variable.
%
% Written by Julian Moosmann.
%
% [im, tif_info] = read_image(filename, filetype, roi, tif_info, shape, dtype, trafo)

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
if nargin < 7
    trafo = '';
end
if nargin < 8
    tifftrafo = 1;
end

%% TODO: CHECK CASES !
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
    case 'tif'
        warning( 'off', 'MATLAB:imagesci:rtifc:missingPhotometricTag' );
        warning( 'off', 'MATLAB:imagesci:tifftagsread:expectedTagDataFormat' );
        warning( 'off', 'imageio:tifftagsread:expectedTagDataFormat' );
        if isempty( tif_info)
            tif_info = imfinfo( filename );
        end
        if ~tifftrafo
            tif_info.Orientation = 0;
        end
        switch tif_info.Orientation
            case 0
                im = imread( filename, 'tif', 'PixelRegion', {[1 tif_info.Height ], [1 tif_info.Width]} );
            case 1
                
                %% MOD: 2020/08/21 for CCD
                if tif_info.Height == tif_info.Width && tif_info.Width == 3056
                    x0 = 1;
                    x1 = tif_info.Height;
                    y0 = 1;
                    y1 = tif_info.Width;
                    if numel( roi ) > 1
                        x0 = max( (tif_info.Height-roi(2)+1), 1 );
                        x1 = min( (tif_info.Height-roi(1)+1), tif_info.Height);                        
                    end
                    im = imread( filename, 'tif', 'PixelRegion', {[y0 y1], [x0 x1 ]} );
                    im = flipud(im);
                else
                    x0 = 1;
                    x1 = tif_info.Height;
                    y0 = 1;
                    y1 = tif_info.Width;
                    
                    %% Implement horizontal ROI
                    %                 if numel( roi ) > 3
                    %                     error( 'Horizontal ROI not yet implemented yet' )
                    %                 end
                    %                 %im = rot90( imread( filename, 'tif', 'PixelRegion', {[x0 x1 ], [y0 y1]} ), 2);
                    %% MOD: 2019-10-04 Ximea camera
                    %                   if numel( roi ) > 1
                    %                     y0 = max( (tif_info.Width-roi(2)+1), 1 );
                    %                     y1 = min( (tif_info.Width-roi(1)+1), tif_info.Width);
                    %                 end
                    %                 if numel( roi ) > 1
                    %                     y0 = max( (tif_info.Width-roi(2)+1), 1 );
                    %                     y1 = min( (tif_info.Width-roi(1)+1), tif_info.Width);
                    %                 end
                    %im = fliplr( rot90( imread( filename, 'tif', 'PixelRegion', {[x0 x1 ], [y0 y1]} ), 1 ) );
                    
                    %% MOD: 2019-12-12
                    if numel( roi ) > 1
                        x0 = max( (tif_info.Height-roi(2)+1), 1 );
                        x1 = min( (tif_info.Height-roi(1)+1), tif_info.Height);
                        
                    end
                    if numel( roi ) == 4
                        y0 = max( roi(3), 1);
                        y1 = min( roi(4), tif_info.Width);
                    end
                    
                    
                    im = rot90( imread( filename, 'tif', 'PixelRegion', {[x0 x1 ], [y0 y1]} ), -1 );
                end
                %im = flipud( rot90( imread( filename, 'tif', 'PixelRegion', {[x0 x1 ], [y0 y1]} ), -1 ) );
                
            case 3
                %% Check vert/hor ROI
                %                 x0 = 1;
                %                 x1 = tif_info.Width;
                %                 y0 = 1;
                %                 y1 = tif_info.Height;
                switch numel( roi )
                    case 0
                        im = flipud( read_tif(filename, tif_info) );
                    case 2
                        im = flipud( read_tif(filename, tif_info, roi ) );
                end
                
            case 4 % Added 2018-11-22, modified 2019-01-30
                %% Check left/right orientation
                x0 = 1;
                x1 = tif_info.Height;
                y0 = 1;
                y1 = tif_info.Width;
                if numel( roi ) > 1
                    x0 = max( (tif_info.Height-roi(2)+1), 1 );
                    x1 = min( (tif_info.Height-roi(1)+1), tif_info.Height);
                end
                if numel( roi ) > 3
                    y0 = max( (tif_info.Width-roi(4)+1), 1 );
                    y1 = min( (tif_info.Width-roi(3)+1), tif_info.Width);
                end
                im = rot90( imread( filename, 'tif', 'PixelRegion', {[x0 x1], [y0 y1]} ), 1);
                
            case 6
                x0 = 1;
                x1 = tif_info.Width;
                y0 = 1;
                y1 = tif_info.Height;
                if numel( roi ) > 1
                    y0 = max( (tif_info.Height-roi(2)+1), 1 );
                    y1 = min( (tif_info.Height-roi(1)+1), tif_info.Height);
                end
                if numel( roi ) > 3
                    x0 = max( (tif_info.Width-roi(4)+1), 1 );
                    x1 = min( (tif_info.Width-roi(3)+1), tif_info.Width);
                end
                im = rot90( imread( filename, 'tif', 'PixelRegion', {[x0 x1], [y0 y1]} ), 2);
                
            case 8 % Added 2018-10-24
                x0 = 1;
                x1 = tif_info.Height;
                y0 = 1;
                y1 = tif_info.Width;
                if numel( roi ) > 1
                    y0 = max( (tif_info.Width-roi(2)+1), 1 );
                    y1 = min( (tif_info.Width-roi(1)+1), tif_info.Width);
                end
                %% Implement horizontal ROI
                if numel( roi ) > 3
                    error( 'Horizontal ROI not yet implemented for TIFF orientation case 4' )
                end
                im = rot90( imread( filename, 'tif', 'PixelRegion', {[x0 x1], [y0 y1]} ), 0);
                
            otherwise
                error( 'TIFF orientation %u not implemented!', tif_info.Orientation )
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'edf'
        im = pmedfread( filename );
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'tif_before20180428'
        if isempty( tif_info)
            tif_info = imfinfo( filename );
        end
        switch numel( roi )
            case 0
                switch tif_info.Orientation
                    case 1
                        % 2018-04-28, FLI writes tif without orientation
                        % flag which defaults to 1
                        im = rot90( fliplr( read_tif(filename, tif_info) ), -1);
                    case 3
                        im = flipud( read_tif(filename, tif_info) );
                end
            case 2
                switch tif_info.Orientation
                    case 1
                        % 2018-04-28, FLI writes tif without orientation flag
                        im = rot90( fliplr( read_tif( filename, tif_info, [(tif_info.Height-roi(2)+1) (tif_info.Height-roi(1)+1)]  ) ) ,-1);
                    case 3
                        im = flipud( read_tif( filename, tif_info, [(tif_info.Height-roi(2)+1) (tif_info.Height-roi(1)+1)] ) );
                end
            otherwise
                error( 'Only vertical ROI supported for fast tif reading.' )
        end
    case {'img', 'dar', 'ref'}
        im = read_dat_jm( filename, roi );
    case 'raw'
        %%rot90 for FLI but not for KIT
        if shape(1) > shape(2)
            im = flipud( read_raw( filename, shape, dtype, roi ) );
            % probably breaks for earlier scans with KIT camera saving raw format, March 2018
            %im = ( read_raw( filename, shape, dtype, roi ) );
        else
            % roi for raw is horizontal roi not verticall!
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

if ~isempty( trafo )
    im = eval( trafo );
end