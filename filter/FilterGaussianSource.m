function filt = FilterGaussianSource( shape_or_image, source_sigma_h_v, dist_source_sample, dist_sample_detector, pixelsize, precision)
% Fourier space filter for multipliation in the detection plane for a
% source with normalized Gaussian profile in real space according to the
% Van Citter-Zernike theorem in the far field. The 2D source size is given
% by source_sigma_h_v corresponding to the RMS standard deviation of the
% Gaussian profile. The Gaussian % source profile is 1/2/pi/sigma^2 * exp(
% - (x/source_sigma_v)^2 -  % (y/source_sigma_h)^2 ). The source is
% situated at a distance % dist_source_sample from the sample and the
% intensity is detected at a % distance dist_sample_detector from the
% sample i.e. at a distance % dist_source_sample+dist_sample_detector from
% the source. 
%
% Written by Julian Moosmann, 2018-02-05. Last version: 2018-02-09

% FT of Gaussian
% F[exp(-a*x^2)](xi) = sqrt(pi/a)*exp(-pi^2/a/xi^2)

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    shape_or_image = [100, 100];
end
if nargin < 2
    % [vert, hor]
    source_sigma_h_v = [6 140] * 1e-6; 
end
if nargin < 3
    dist_source_sample = 82.7; % m
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
if numel( source_sigma_h_v ) == 1
    source_sigma_h_v = [source_sigma_h_v, source_sigma_h_v];
end

% source size
source_sigma_h = source_sigma_h_v(1);
source_sigma_v = source_sigma_h_v(2);

% 1D grid
% Check factor 1/2
x = 1 / 2 / pixelsize * FrequencyVector(shape_or_image(2), precision, 1);
y = 1 / 2 / pixelsize *  FrequencyVector(shape_or_image(1), precision, 1);

% 2D grid
[xx, yy] = meshgrid( x, y);

% Fourier space filter
%filt = exp( - ( pi^2 * dist_sample_detector ) / ( 4 * log(2) * dist_source_sample) * ( (source_sigma_h * xx).^2 + (source_sigma_v * yy).^2 ) );
filt = exp( - 2 * ( pi * dist_sample_detector  / dist_source_sample )^2 * ( (source_sigma_h * xx).^2 + (source_sigma_v * yy).^2 ) );

% %% Blurring
% 
% % source: standard deviation via FWHM for a Gaussian profile
% source_sigma_v = source_sigma_h / ( 2 * sqrt( 2 * log(2) ) );
% source_sigma_h = source_sigma_v / ( 2 * sqrt( 2 * log(2) ) );
% 
% % source blurring
% blur_sigma_v = dist_sample_detector / dist_source_sample * source_sigma_v;
% blur_sigma_h = dist_sample_detector / dist_source_sample * source_sigma_h;
% blur_sigma_h_v = [blur_sigma_h, blur_sigma_v];
% 
% % "cut off" frequency
% blur_cutoff_frequ_h_v = 1 / sqrt(2) / pi ./ [blur_sigma_h, blur_sigma_v];
% 
% % Fourier space filter
% filt_b =  exp( - 2 * pi^2 *  dist_source_sample / dist_sample_detector *( (blur_sigma_v*xx).^2 + (blur_sigma_h*yy).^2 ) );
% 
% 
