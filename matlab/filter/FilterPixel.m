function [im_int, threshold_hot, threshold_dark] = FilterPixel( im_int, par, gpu_index )
% Filter hot, dark, or dead pixels of the 2D input image. Filtered pixels
% are substituted by the median values of their neighboorhood. Hot and dark
% pixels are identified using the ratio of the image and the median
% filtered images which has to be bigger than a certain threshold.
%
% im : 2D array, image to be filtered
% par : parameter struct with fields:
%   threshold_hot : positive scalar, default: 0.01. Values < 1 are
%       interpreted as relative numbers of pixels to be filtered. Values >=
%       1 are interpreted as thresholds. Pixels of the matrix
%       'im./medfilt2( im )' with a threshold.
%   threshold_dark : positive scalar, default: 0. Values < 0.5 are
%       interpreted as relative numbers of pixels to be filtered. Values >=
%       0.5 are interpreted as thresholds. Pixels of the matrix
%       'im ./ medfilt2( im )' with a threshold.
%   medfilt_neighboorhood: 2-vector, default: [3 3]. Neighborhood of the
%       median filter applied to the image.
%   filter_dead_pixel : logical, default: true.
%   filter_Inf : logical, default: true.
%   filter_NaN : logical, default: true.
%   verbose : logical, default: 0.
%   use_gpue :  logical, default
% gpu_index : 1D-vector of integer indices of the GPU devices to be used.
% indices start from 1
%
% If threshold_hot, threshold_dark, and  filter_dead_pixel are all zero,
% the input image is returned.
%
% Calculating the threshold from the prescription to filter X % of all pixel
% is more compuationally extensive compared to the case where the values are
% given directly, and thus should be avoided for a large amount of data.
%
% Note: Printing information about the unfiltered and filtered image
% produces additonal workload and should be turned off for speed.
%
% Written by Julian Moosmann.
%
% [im, thresholds, percentages] = FilterPixelGPU( im, par )

%warning( 'off', 'parallel:gpu:device:DeviceDeprecated' )

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    par.use_gpu = 0;
    par.verbose = 0;
end
if nargin < 3
    gpu_index = [];
end
use_gpu = par.use_gpu;

if use_gpu && ~isempty( gpu_index )
    p = parallel.gpu.GPUDevice.current();
    cur_gpu_index = p.Index;
    if cur_gpu_index ~= gpu_index
        gpuDevice( gpu_index );
    end
end
verbose = assign_from_struct( par, 'verbose', 0 );
if verbose
    t = toc;
    im_min  = min( im_int(:) );
    im_max  = max( im_int(:) );
    im_mean = mean( im_int(:) );
    im_std  = std( single( im_int(:) ) );
    inp_class = class( im_int );
    if use_gpu
        gpu = parallel.gpu.GPUDevice.current;
        mem0 = gpu.AvailableMemory;
    end
end
threshold_hot = assign_from_struct( par, 'threshold_hot', 0.01 );
threshold_dark = assign_from_struct( par, 'threshold_dark', 0, 1 );
medfilt_neighboorhood = assign_from_struct( par, 'medfilt_neighboorhood', [3 3] );
filter_dead_pixel = assign_from_struct( par, 'filter_dead_pixel', 1 );
filter_Inf = assign_from_struct( par, 'filter_Inf', 1 );
filter_NaN = assign_from_struct( par, 'filter_NaN', 1 );

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%gpu = parallel.gpu.GPUDevice.current;
%fprintf( ' %u', gpu.Index );

mfnh1 = medfilt_neighboorhood(1);
mfnh2 = medfilt_neighboorhood(2);

% Return original image if threshhold are zero.
if isequal( threshold_hot, 0) && isequal( threshold_dark, 0) && isequal( filter_dead_pixel, 0 )
    return;
end

% Create GPU array
if use_gpu
    im = gpuArray( im_int );
else
    im = im_int;
end

% Pad array since medfilt on GPU has no padding option
im = padarray( im, [mfnh1 mfnh2], 'symmetric', 'both' );

% Info
num_pix = numel( im );
num_dead = 0;
num_hot  = 0;
num_dark = 0;
num_NaN = 0;
num_Inf = 0;

% Filter mask
if use_gpu
    mask = zeros( size(im), 'logical','gpuArray' );
else
    mask = zeros( size(im), 'logical' );
end

% Detect INFs
if isfloat( im_int ) && filter_Inf
    mask = mask | isinf( im );
    num_Inf = sum( mask(:) );
end

% Detect NANs
if isfloat( im_int ) && filter_NaN
    mask = mask | isnan( im );
    num_NaN = sum( mask(:) ) - num_Inf;
end

% Dead pixel mask
if filter_dead_pixel > 0
    mask = mask | im <= 0;
    num_dead = sum( mask(:) ) - num_Inf - num_NaN;
    % Return if dark/hot pixels are not filtered & no dead pixels detected
    if threshold_hot == 0 && threshold_dark == 0 && num_dead == 0
        return
    end
end

% Replace Infs, NaNs, and dead pixel before median filtering
if filter_Inf || filter_NaN || filter_dead_pixel
    im_mean2 = mean2( im(~mask) );
    im(mask) = im_mean2;
end

% Median filtered image
im_med = medfilt2( im, medfilt_neighboorhood );

% Ratio of image and median filtered image
R = single( im );
R = R ./ single( im_med );
%% remove R_sorted to reduce memory
if (threshold_hot > 0 && threshold_hot < 1) || (threshold_dark > 0 && threshold_dark < 0.5)
    R_sorted = sort( R(:) );
end

% Hot pixel mask
if threshold_hot > 0
    % Determine threshold_hot if not given explixitly as a value >= 1par
    if threshold_hot < 1
        x = floor( num_pix * ( 1 - threshold_hot ) );
        threshold_hot = gather( R_sorted(x) );
    end
    % Combine hot pixel mask and previos mask
    mask = mask | R > threshold_hot;
    num_hot = sum( mask(:) ) - num_dead;
end

% Dark pixel mask
if threshold_dark > 0
    % Determine threshold_dark if not given explixitly as a value >= 0.5
    if threshold_dark < 0.5
        x = floor( num_pix * threshold_dark );
        threshold_dark = gather( R_sorted(x) );
    end
    mask = mask | R < threshold_dark;
    num_dark = sum( mask(:) ) - num_dead - num_hot;
end

% Replace pixels to be filtered
im(mask) = im_med(mask);

% Retrieve form GPU
if use_gpu
    im_int = gather( im(mfnh1+1:end-mfnh1,mfnh2+1:end-mfnh2) );     
else
    im_int = im(mfnh1+1:end-mfnh1,mfnh2+1:end-mfnh2) ;
end

%% Print info
if verbose
    t = toc - t ;
    im_filt_min  = min( im_int(:) );
    im_filt_max  = max( im_int(:) );
    im_filt_mean = mean( im_int(:) );
    im_filt_std  = std( single( im_int(:) ) );
    fprintf( '\n class : %s (Note: Inf/NaN detection works only for float arrays)', inp_class );
    fprintf( '\n shape padded : %u %u', size( im ) )
    fprintf( '\n number of pixels: %u', num_pix );
    if filter_dead_pixel > 0
        fprintf( '\n number of dead pixels: %9u (%3f %%)', num_dead, 100*num_dead/num_pix )
    end
    if threshold_hot > 0
        fprintf( '\n number  of hot pixels: %9u (%3f%%), threshold: %9g', num_hot, 100 * num_hot/num_pix, threshold_hot );
    end
    if threshold_dark > 0
        fprintf( '\n number of dark pixels: %9u (%3f%%), threshold: %9g', num_dark, 100 * num_dark/num_pix, threshold_dark );
    end
    if num_Inf > 0
        fprintf( '\n number of Infs: %9u (%3f%%)', num_Inf, 100 * num_Inf / num_pix );
    end
    if num_NaN > 0
        fprintf( '\n number of NANs: %9u (%3f%%)', num_NaN, 100 * num_NaN / num_pix );
    end
    fprintf( '\n before filter: [Min Max Mean Std] = [%9g %9g %9g %9g]', im_min, im_max, im_mean, im_std );
    fprintf( '\n after filter : [Min Max Mean Std] = [%9g %9g %9g %9g]', im_filt_min, im_filt_max, im_filt_mean, im_filt_std );
    if use_gpu
        mem1 = gpu.AvailableMemory;
        memt = gpu.TotalMemory;
        
        fprintf( '\n gpu device index : %u', gpu.Index )
        fprintf( '\n gpu memory at start: %f (%g)', mem0 / 1024^2, mem0 / memt )
        fprintf( '\n gpu memory at end  : %f (%g)', mem1 / 1024^2, mem1 / memt )
        fprintf( '\n gpu memory diff.   : %f (%g)', (mem0 - mem1) / 1024^2, (mem0 - mem1) / memt )
    end
    fprintf( '\n Elapsed time : %g s', t  )
    fprintf( '\n' )
end

clear im im_med mask R
if exist('R_sorted','var')
    clear R_sorted
end
if exist('im_mean2','var')
    clear im_mean2
end
if use_gpu
    clear num_dead num_hot num_NaN num_dark num_pix
end

