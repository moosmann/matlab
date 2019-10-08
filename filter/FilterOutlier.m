function [im, im_low, im_high] = FilterOutlier( im, threshold, bit_conversion, return_threshs_only, verbose)
% Remove outliers (and optionally convert to unsigned integer.
%
% im : 2D/3D-array
% threshold : 1- or 2-component vector. Outliers to be filtered. If scalar
%   lowest ands highest outliers are filtered alike.
% bit_conversion : str, default: ''. if not empty convert to unsigned 8-,
%   16-, 32-, or 64-bit integer. if empty: rescale, but do not convert.
% return_threshs_only : bool, default:0. If only thresholds are desired.
% verbose : bool, default:1. Print information.
%
% Written by J. Moosmann.
%
% [im, im_low, im_high] = FilterOutlier( im, threshold, bit_conversion, return_threshs_only, verbose)

%% Default arguments
if nargin < 2
    threshold = 0.03;
end
if nargin < 3
    bit_conversion = '';
end
if nargin < 4
    return_threshs_only = 0;
end
if nargin < 5
    verbose = 0;
end

%% Thresholds
if numel( threshold ) == 1
    threshold = [threshold, threshold];
end

% Sort image values
im_sorted = sort( im(:) );
PrintVerbose( verbose, '\n image : min, max = %f, %f', im_sorted(1), im_sorted(end) )
% Lower threshold and index
ind_low = ceil( threshold(1) * numel( im_sorted ) );
im_low = im_sorted( ind_low );
% Higher threshold and index
ind_high = floor( threshold(2) * numel( im_sorted ) );
im_high = im_sorted(end-ind_high);
PrintVerbose( verbose, '\n threshold low %g : %f', threshold(1), im_low )
PrintVerbose( verbose, '\n threshold high %g : %f', threshold(2), im_high )

if ~return_threshs_only
    
    %% Filter
    if nargout > 0
        im( im < im_low) = im_low;
        im( im > im_high) = im_high;
        PrintVerbose( verbose, '\n Cropping dynamic range.' )
    end
    
    %% bit conversion
    if ~isempty( bit_conversion )
        switch bit_conversion
            case {'8bit', 8, 'uint8' }
                im =  uint8( (2^8 - 1) * (im - im_low) / (im_high - im_low ) );
            case {'16bit', 16, 'uint16' }
                im =  uint16( (2^16 - 1) * (im - im_low) / (im_high - im_low ) );
            case {'32bit', 32, 'uint32'}
                im =  uint32( (2^32 - 1) * (im - im_low) / (im_high - im_low ) );
            case {'64bit', 64, 'uint64'}
                im =  uint64( (2^64 - 1) * (im - im_low) / (im_high - im_low ) );
        end
        PrintVerbose( verbose, '\n Bit conversion.' )
    else
        im = (im - im_low) / (im_high - im_low );
    end
end
PrintVerbose( verbose, '\n' )
