function [im, tif_info] = read_image( filename, par, no_raw_roi )
% Read 'tif', 'img', 'dar', 'ref', and 'edf' images and all formats that
% are MATLAB's imread function can read without additional arguments than file
% format.
%
% ARGUMENT:
% % filename: string. Absolute filename
% par: parameter struct with fields as
%
% FIELDS of 'par':
% im_format: string. Suffix of image file format. Default: ''
% raw_roi : Region of interest
% tif_info : image info struct
% im_shape_raw : required for raw images
% im_trafo : string, default ''. string which will be evaluated on the read
%   image, e.g. 'rot90(im)'. Note that 'im' is the mandatory variable.
% tifftrafo : use tiff orientagion tag to determine the image
%   transformation.
%
% Written by Julian Moosmann.
%
% [im, tif_info] = read_image( filename, par )

% old: function [im, tif_info] = read_image(filename, dtype, raw_roi, tif_info, im_shape_raw, dtype, im_trafo, tifftrafo)

%% Default arguments
if nargin < 2
    par = struct;
end
if nargin < 3
    no_raw_roi = 0;
end
assign_from_struct( par, 'im_format',  '');
assign_from_struct( par, 'raw_roi', [] );
assign_from_struct( par, 'tif_info', [] );
assign_from_struct( par, 'im_shape_raw', [] );
assign_from_struct( par, 'dtype',  '');
assign_from_struct( par, 'im_trafo', '' );
assign_from_struct( par, 'tifftrafo',  1);


%% TODO: CHECK CASES !
%% TODO: raw_roi support for edf and tiff files
%% TODO: check flip/rotation/transponation for raw data

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty( im_format )
    [~, ~, im_format] = fileparts(filename);
end
if strcmp( im_format(1), '.')
    im_format(1) = [];
end
if no_raw_roi
    raw_roi = [];
end

switch lower( im_format )
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
                    if numel( raw_roi ) > 1
                        x0 = max( (tif_info.Height-raw_roi(2)+1), 1 );
                        x1 = min( (tif_info.Height-raw_roi(1)+1), tif_info.Height);                        
                    end
                    im = imread( filename, 'tif', 'PixelRegion', {[y0 y1], [x0 x1 ]} );
                    im = flipud(im);
                else
                    % hor/vert ROI works for Ximea 50 MP, tested
                    % 2021/11/16, tested with imsc1
                    x0 = 1; % top
                    x1 = tif_info.Height; % bottom
                    y0 = 1; % left
                    y1 = tif_info.Width; % right
                    % MOD: 2019-12-12
                    if numel( raw_roi ) > 1
                        x0 = max( (tif_info.Height-raw_roi(2)+1), 1 );
                        x1 = min( (tif_info.Height-raw_roi(1)+1), tif_info.Height);
                    end
                    if numel( raw_roi ) == 4
                        y0 = max( raw_roi(3), 1);
                        y1 = min( raw_roi(4), tif_info.Width);
                    end
                    im = rot90( imread( filename, 'tif', 'PixelRegion', {[x0 x1 ], [y0 y1]} ), -1 );
                end
            case 3
                % READ TIFF about twice as fast as imread
                %                 switch numel( raw_roi )
                %                     case 0
                %                         im = flipud( read_tif(filename, tif_info) );
                %                     case 2
                %                         im = flipud( read_tif(filename, tif_info, raw_roi ) );
                %                 end
                
                % Check for KIT camera
                if tif_info.Width == 5120
                    x0 = 1;
                    x1 = tif_info.Height;                    
                    y0 = 1;
                    y1 = tif_info.Width;                    
                    if numel( raw_roi ) > 1
                        x0 = max( (raw_roi(1)+1), 1 );
                        x1 = min( (raw_roi(2)+1), tif_info.Height);
                    end
                else
                    %% MOD: 2020-09-03 Ximea sCMOS
                    x0 = 1;
                    x1 = tif_info.Height;
                    %                 x1 = tif_info.Width;
                    y0 = 1;
                    y1 = tif_info.Width;
                    %                 y1 = tif_info.Height;
                    if numel( raw_roi ) > 1
                        x0 = max( (tif_info.Height-raw_roi(2)+1), 1 );
                        x1 = min( (tif_info.Height-raw_roi(1)+1), tif_info.Height);
                    end
                    if numel( raw_roi ) == 4
                        y0 = max( raw_roi(3), 1);
                        y1 = min( raw_roi(4), tif_info.Width);
                    end
                end
                try
                    im = rot90( imread( filename, 'tif', 'PixelRegion', {[x0 x1 ], [y0 y1]} ), 1 );
                catch
                    im = zeros([y1 - y0 + 1, x1 - x0 + 1], 'uint16');
                end
                
            case 4 % Added 2018-11-22, modified 2019-01-30
                %% Check left/right orientation
                x0 = 1;
                x1 = tif_info.Height;
                y0 = 1;
                y1 = tif_info.Width;
                if numel( raw_roi ) > 1
                    x0 = max( (tif_info.Height-raw_roi(2)+1), 1 );
                    x1 = min( (tif_info.Height-raw_roi(1)+1), tif_info.Height);
                end
                if numel( raw_roi ) > 3
                    y0 = max( (tif_info.Width-raw_roi(4)+1), 1 );
                    y1 = min( (tif_info.Width-raw_roi(3)+1), tif_info.Width);
                end
                im = rot90( imread( filename, 'tif', 'PixelRegion', {[x0 x1], [y0 y1]} ), 1);
                
            case 6
                x0 = 1;
                x1 = tif_info.Width;
                y0 = 1;
                y1 = tif_info.Height;
                if numel( raw_roi ) > 1
                    y0 = max( (tif_info.Height-raw_roi(2)+1), 1 );
                    y1 = min( (tif_info.Height-raw_roi(1)+1), tif_info.Height);
                end
                if numel( raw_roi ) > 3
                    x0 = max( (tif_info.Width-raw_roi(4)+1), 1 );
                    x1 = min( (tif_info.Width-raw_roi(3)+1), tif_info.Width);
                end
                im = rot90( imread( filename, 'tif', 'PixelRegion', {[x0 x1], [y0 y1]} ), 2);
                
            case 8 % Added 2018-10-24
                x0 = 1;
                x1 = tif_info.Height;
                y0 = 1;
                y1 = tif_info.Width;
                if numel( raw_roi ) > 1
                    y0 = max( (tif_info.Width-raw_roi(2)+1), 1 );
                    y1 = min( (tif_info.Width-raw_roi(1)+1), tif_info.Width);
                end
                %% Implement horizontal ROI
                if numel( raw_roi ) > 3
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
        switch numel( raw_roi )
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
                        im = rot90( fliplr( read_tif( filename, tif_info, [(tif_info.Height-raw_roi(2)+1) (tif_info.Height-raw_roi(1)+1)]  ) ) ,-1);
                    case 3
                        im = flipud( read_tif( filename, tif_info, [(tif_info.Height-raw_roi(2)+1) (tif_info.Height-raw_roi(1)+1)] ) );
                end
            otherwise
                error( 'Only vertical ROI supported for fast tif reading.' )
        end
    case {'img', 'dar', 'ref'}
        im = read_dat_jm( filename, raw_roi );
    case 'raw'
        %%rot90 for FLI but not for KIT
        if im_shape_raw(1) > im_shape_raw(2)
            im = flipud( read_raw( filename, im_shape_raw, dtype, raw_roi ) );
            % probably breaks suport for previou scans with KIT camera saving raw format, March 2018
            %im = ( read_raw( filename, im_shape_raw, dtype, raw_roi ) );
        else
            % roi for raw is horizontal roi not verticall!
            if isempty( raw_roi )
                im = flipud(rot90(read_raw( filename, im_shape_raw, dtype, [] ), -1));
            else
                im = read_raw( filename, im_shape_raw, dtype, [] );
                im = im(raw_roi(1):raw_roi(2),:);
                im = flipud(rot90(im, -1));
            end
        end
    otherwise
        im = imread( filename, im_format )';
end

if ~isempty( im_trafo )
    im = eval( im_trafo );
end
