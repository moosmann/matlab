function [filt, blur_sigma_h_v, blur_cutoff_frequ_h_v] = FilterGaussianSource( shape_or_image, source_size_h_v, dist_source_sample, dist_sample_detector, pixelsize, precision)
% Gaussian filter in Fourier space for a normalized Gaussian source profile in
% real space with a source size (s_x,s_y) at FWHM and a profile
% 1/2/pi/sigma^2 * exp( - (x/blur_sigma_v)^2 - (y/blur_sigma_h)^2 ), situated at a
% distance R from the sample and detected at a distance z from the sample
% ie a distance z + R from the source.
%
% Written by Julian Moosmann, 2018-02-05

if nargin < 1
    shape_or_image = [100, 100];
end
if nargin < 2
    % [vert, hor]
    source_size_h_v = [6 140] * 1e-6; 
end
if nargin < 3
    dist_source_sample = 90; % m
end
if nargin < 4
    dist_sample_detector = 1; % m
end
if nargin < 5
    pixelsize = 1e-6; % m
end
if nargin < 6
    precision = 'double';
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if numel( shape_or_image ) > 2
    shape_or_image = size( shape_or_image );
end
if numel( source_size_h_v ) == 1
    source_size_h_v = [source_size_h_v, source_size_h_v];
end

% source: standard deviation via FWHM for a Gaussian profile
source_sigma_v = source_size_h_v(1) / ( 2 * sqrt( 2 * log(2) ) );
source_sigma_h = source_size_h_v(2) / ( 2 * sqrt( 2 * log(2) ) );

% source blurring
blur_sigma_v = dist_sample_detector / dist_source_sample * source_sigma_v;
blur_sigma_h = dist_sample_detector / dist_source_sample * source_sigma_h;
blur_sigma_h_v = [blur_sigma_h, blur_sigma_v];

% "cut off" frequency
blur_cutoff_frequ_h_v = 1 / sqrt(2) / pi ./ [blur_sigma_h, blur_sigma_v];

% 1D grid
x = 1 / 2 / pixelsize * FrequencyVector(shape_or_image(2), precision, 1);
y = 1 / 2 / pixelsize *  FrequencyVector(shape_or_image(1), precision, 1);
% 2D grid
[xx, yy] = meshgrid( x, y);

% Fourier space filter (prefactor cancels because of normalization)
filt =  exp( - 2 * pi^2 *  dist_source_sample / dist_sample_detector *( (blur_sigma_v*xx).^2 + (blur_sigma_h*yy).^2 ) );
