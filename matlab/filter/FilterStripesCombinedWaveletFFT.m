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
% wname : string, wavelet type. default: 'db25'. 'bior2.8'
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

% Daubechies
% 	'db1' or 'haar', 'db2', ... ,'db10', ... , 'db45'
% 
% Coiflets
% 	'coif1', ... , 'coif5'
% 
% Symlets
% 	'sym2', ... , 'sym8', ... ,'sym45'
% 
% Fejer-Korovkin filters
% 	'fk4', 'fk6', 'fk8', 'fk14', 'fk22'
% 
% Discrete Meyer
% 	'dmey'
% 
% Biorthogonal
% 	'bior1.1', 'bior1.3', 'bior1.5'
% 'bior2.2', 'bior2.4', 'bior2.6', 'bior2.8'
% 'bior3.1', 'bior3.3', 'bior3.5', 'bior3.7'
% 'bior3.9', 'bior4.4', 'bior5.5', 'bior6.8'
% 
% Reverse Biorthogonal
% 	'rbio1.1', 'rbio1.3', 'rbio1.5'
% 'rbio2.2', 'rbio2.4', 'rbio2.6', 'rbio2.8'
% 'rbio3.1', 'rbio3.3', 'rbio3.5', 'rbio3.7'
% 'rbio3.9', 'rbio4.4', 'rbio5.5', 'rbio6.8'

%% TODO: add direction option

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    im = zeros(100,100);
    im(10:4:end-9) = 1;
end
if nargin < 2
    dec_levels = 1:7;
end
if nargin < 3
   wname = 'db25'; 
end
if nargin < 4
    sigma = 2.4;
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[d1, d2] = size( im );
dec_num = max( dec_levels );

% wavelet decomposition
Ch{dec_num} = [];
Cv{dec_num} = [];
Cd{dec_num} = [];
for ii = 1:dec_num
    [im, Ch{ii}, Cv{ii}, Cd{ii}] = dwt2( im, wname);    
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

im_filt = im_filt(1:d1,1:d2);

return
