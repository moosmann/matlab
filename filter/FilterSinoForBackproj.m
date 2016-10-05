function sino = FilterSinoForBackproj(sino, direction, filt_type, pad_method, filt_len, filt_frequ_cutoff)
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

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    direction = 1;
end
if nargin < 3
    filt_type = 'Ram-Lak';
end
if nargin < 4
    pad_method = 0;
end
if nargin < 5
    filt_len = 'twice';
end
if nargin < 6
    filt_frequ_cutoff = 1;
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
            error(sprintf('Length of filter FILT_LEN must not be smaller than %u', len))
        end        
    end

%     % before
%     pad_vec = zeros([1, ndims(sino)] );
%     pre_size = floor( len / 2 );
%     pad_vec(direction) = pre_size;
%     sino = padarray(sino, pad_vec, pad_method , 'pre');
%     % after
%     pad_vec = zeros([1, ndims(sino)] );
%     post_size = ceil( len / 2 );
%     pad_vec(direction) = post_size;
%     sino = padarray(sino, pad_vec,  pad_method, 'post');
    
    pad_vec = zeros([1, ndims(sino)] );
    pad_vec(direction) = filt_len - len;
    sino = padarray(sino, pad_vec,  pad_method, 'post');
end
disp(size(sino))
% Design filter
filt = iradonDesignFilter(filt_type, filt_len, filt_frequ_cutoff);
filt_vec = ones([1, ndims(sino)] );
filt_vec(direction) = filt_len;
filt = reshape( filt, filt_vec);

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
%             sino = sino(pre_size + (1:len) , :, : );
            sino( len+1:end, :, :) = [];
        case 2
%             sino = sino(:, pre_size + (1:len), : );
            sino( :, len+1:end, :) = [];
        case 3
%             sino = sino(:, :, pre_size + (1:len) );
            sino( :, :, len+1:end) = [];
    end

end
