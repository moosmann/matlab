function [im_filt] = FilterStripesCombinedWaveletFFT(im, dec_levels, wname, sigma)
% Combined wavelet and Fourier analysis for the elimination or vertical
% stripes artifacts.
%
% This algorithm employs three
% parameters for the filtering process: the highest decomposition level L
% (’dec_levels’), the wavelet type (’wname’) and the damping factor s (’sigma’)
% from Eq.9. In our approach, s is kept uniform for each decomposition
% level l ∈ {0, · · · ,L}. Alternatively, s could also be adjusted
% according to the resolution at each l.
%
% PARAMETER:
% im : 2D input image to be stripe filtered
% dec_levels : 1D-vector of integers in [1, max_level]. decompositon level 
% wname : string, wavelet type. default: 'db25'
% sigma : scalar or 1D-vector with numel(sigma)=max(dec_levels), damping
% factor, if scalar same sigma is used for all levels.
%
% Taken from Münch, B.; Trtik, P.; Marone, F. & Stampanoni, M. Stripe and
% ring artifact removal with combined wavelet - Fourier filtering Opt.
% Express, OSA, 2009, 17, 8567-8591.
% 
% Modified by Julian Moosmann, 2017-05-30.
%
% [im_filt] = FilterStripesCombinedWaveletFFT(im, dec_levels, wname, sigma)

%% TODO: add direction option

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 3
    im = zeros(100,100);
    im(10:4:end-9) = 1;
end
if nargin < 2
    dec_levels = 2:5;
end
if nargin < 3
   wname = 'db25'; 
end
if nargin < 4
    sigma = 2.4;
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dec_num = max( dec_levels );

% wavelet decomposition
Ch{dec_num} = [];
Cv{dec_num} = [];
Cd{dec_num} = [];
for ii = 1:dec_num
    [im, Ch{ii}, Cv{ii}, Cd{ii}] = dwt2(im, wname);
end

% FFT transform of horizontal frequency bands
if numel( sigma ) == 1
    sigma = sigma * ones( dec_num, 1);
end
for ii = dec_levels
    
    % FFT
    fCv = fftshift( fft( Cv{ii} ) );
    [my, mx] = size( fCv );
    
    % damping of vertical stripe information
    damp = 1 - exp( -(-floor(my/2):-floor(my/2)+my-1).^2 / (2 * sigma(ii)^2) );
    fCv = fCv .* repmat( damp', 1, mx );
    
    % inverse FFT
    Cv{ii} = ifft( ifftshift( fCv ) );
end

% wavelet reconstruction
im_filt = im;
for ii = dec_num:-1:1
    im_filt = im_filt(1:size(Ch{ii},1),1:size(Ch{ii},2));
    im_filt = idwt2( im_filt, Ch{ii}, Cv{ii}, Cd{ii}, wname);
end

return
