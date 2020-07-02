function sino = FilterSinoForBackproj(sino, direction, filt_type, pad_method, filt_len, filt_frequ_cutoff, butterworth)
% Filter designed to apply to sinogram before backprojection.
%
% sino: 2-or-3-dimensional array
% direction: scalar. Default: 1. Dimesnion along the filter should be
% applied
% filt_type: string. Default: 'Ram-Lak'. Options: 'none', 'Ram-Lak',
% 'Shepp-Logan', 'Cosine', 'Hamming', 'Hann'. See iradonDesingFilter for
% Details.
% pad_method: scalar, 'none', 'replicate', or 'symmetric'. Default: 0.
%
% Currently only even number of pixels in the filtering direction is
% possible. iradonDesignFilter does not yet support odd pixels.
%
% Written by Julian Moosmann
% First version: 2016-09-29. Last modification: 2016-09-30
%
% sino = FilterSinoForBackproj(sino, direction, filt_type)

%% TODO: improve documentantion
%% TODO: improve how padding parameters are provided

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    direction = 1;
end
if nargin < 3
    filt_type = 'Ram-Lak';
end
if nargin < 4
    pad_method = 'symmetric';
end
if nargin < 5
    filt_len = 'twice';
end
if nargin < 6
    filt_frequ_cutoff = 1;
end
if nargin < 7
    butterworth = 0.5;
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if direction == 0
    return
end

len = size( sino, direction); 

% Padding method
if strcmpi(pad_method, 'none')
    filt_len = len;
else
    % Padding size
    if strcmpi(filt_len, 'twice')
        filt_len = 2 * len;
    elseif strcmpi(filt_len, 'nextpow2')
        filt_len = max(64,2^nextpow2(2*len));
    else
        if filt_len < len
            error('Length of filter FILT_LEN must not be smaller than %u', len)
        end        
    end
    
    pad_vec = zeros([1, ndims(sino)] );
    pad_vec(direction) = filt_len - len;
    sino = padarray(sino, pad_vec,  pad_method, 'post');
end

% Design filter
filt = iradonDesignFilter(filt_type, filt_len, filt_frequ_cutoff);
filt_vec = ones([1, ndims(sino)] );
filt_vec(direction) = filt_len;
filt = reshape( filt, filt_vec);

% Butterworth filter
if butterworth > 0
    [b, a] = butter(6, butterworth);
    bw = freqz(b, a, numel(filt) );
    filt = filt .* bw;
end

% Fourier transform
sino = fft( sino, [], direction); % zero frequency is at first position 

% Filter multiplication
sino = bsxfun(@times, sino, filt);
    
% Inverse Fourier transform
sino = real( ifft( sino , [], direction, 'symmetric') );

% Crop to original size
if ~strcmpi(pad_method, 'none')
    switch direction
        case 1
            sino( len+1:end, :, :) = [];
        case 2
            sino( :, len+1:end, :) = [];
        case 3
            sino( :, :, len+1:end) = [];
    end

end
