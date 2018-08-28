function vol = read_multitif( filename, images_to_read, PixelRegion, tiff_info, verbose )
% Read multi-tiff file.
%
% ARGUMENTS
% filename : string
% images_to_read: vector of integer indices to read. Default: [] read all.
% PixelRegion : cell of 2- or 3- component vector as for MATLAB's imread
%   function. Default: {} reads full region.
% tiff_info : struct from MATLAB's imfinfo. Default: [] reads from file.
% verbose : bool, print info. Default: 1
%
% Written by J. Moosmann, last version: 2018-08-24
%
% vol = Readmultitif( filename, images_to_read, PixelRegion, tiff_info, verbose )

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    filename = 'radio.tif';
end
if nargin < 2
    images_to_read = [];
end
if nargin < 3
    PixelRegion = {};
end
if nargin < 4
    tiff_info = [];
end
if nargin < 5
    verbose = 1;
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%warning('off','MATLAB:tifftagsread:tiffTag:wrongTagDataFormat')

t = toc;

PrintVerbose( verbose, 'Read multi-tif file' );

% Info
if isempty( tiff_info )
    tiff_info = imfinfo(filename,'tif');
end
num_tiff = numel(tiff_info);
switch tiff_info(1).BitDepth
    case 32
        dtype = 'single';
    case 64
        dtype = 'double';
end

% Tiffs to read
if isempty( images_to_read )
    images_to_read = 1:num_tiff;
end

num_read = numel(images_to_read);

% ROI defined by PixelRegion
if isempty( PixelRegion )
    shape = [tiff_info(1).Height, tiff_info(1).Width, num_read];
    PixelRegion = { [1 shape(1) 1], [1 shape(2) 1]};
else
    switch numel( PixelRegion{1} )
        case 2
            num_rows = PixelRegion{1}(2) - PixelRegion{1}(1) + 1;
        case 3
            num_rows = numel( PixelRegion{1}(1):PixelRegion{1}(2):PixelRegion{1}(3) );
    end
    switch numel( PixelRegion{2} )
        case 2
            num_cols = PixelRegion{2}(2) - PixelRegion{2}(1) + 1;
        case 3
            num_cols = numel( PixelRegion{2}(1):PixelRegion{2}(2):PixelRegion{2}(3) );
    end
    shape = [num_rows, num_cols, num_read];
end
PrintVerbose( verbose, ', vol shape: [%u %u %u]', shape )

% Preallocate memory
vol = zeros( shape, dtype );

% Read
for nn = 1:num_read
    vol(:,:,nn) = imread( filename, 'Index', images_to_read(nn), 'Info', tiff_info, 'PixelRegion', PixelRegion);
end

PrintVerbose( verbose, ', done in %.1f s = %.3f min\n', toc -t, (toc - t)/60)
