function [mask, half_mask, ctf_filter, quasiparticle_filter] = BinaryMaskCTF(filter_size, energy, distance, pixelsize, regularization_parameter, binary_filter_threshold, phase_shift)
% Create a binary where bands centered at the zero crossings of the
% contrast-transfer function (CTF) are set to zero.
%
% ARGUMENTS
% see function: 'PhaseFilter'
%
% Written by Julian Moosmann, First verion: 2017-02-08. Last: 2017-02-08

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%regularization_parameter = 2.5;
%binary_filter_threshold = 0.2;
precision = 'single';
if nargin < 7
    phase_shift = 0;
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda    = 6.62606896e-34*299792458/(energy/1000*1.60217733e-16);
ArgPrefac = pi*lambda*distance/pixelsize^2;

%% Fourier coordinates
% 1D
xi  = FrequencyVector(filter_size(2),precision,1);
eta = FrequencyVector(filter_size(1),precision,1);
% 2D
[sinArg, sinxiquad]   = meshgrid(xi,eta);
sinArg = ArgPrefac*(sinArg.^2 + sinxiquad.^2);
sinxiquad = sin(sinArg - phase_shift);

%% CTF filter
ctf_filter = 1/2*sign(sinxiquad)./(abs(sinxiquad)+10^-regularization_parameter);

%% Quasiparticle filter
mask = sinArg > pi/2  &  abs(sinxiquad) < binary_filter_threshold;
quasiparticle_filter = ctf_filter;
quasiparticle_filter( mask ) = 0;

%% Binary mask
mask = 1 - fftshift( mask );

%% Half binary mask
half_mask = mask;
half_mask(:, ceil(size(mask,2)/2):end) = 1;