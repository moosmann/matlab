function [int] = propagate(phase_shift, absorption, EnergyDistancePixelsize, padding, padding_method, filt, verbose)
% Intensity according Fresnel propagation for a monochromatic plane wave of
% energy 'EnergyDistancePixelsize'(1), a sample to detector distance
% 'EnergyDistancePixelsize'(2), and detector pixelsize of
% 'EnergyDistancePixelsize'(3) using the 2D maps 'phase_shift'
% and 'absorption' as input. Here 'absorption' is two times the linear
% attenuation, such that the input wave I(z=0) at zero distance is given by
% I(z=0) = exp( i*'phase_shift' - 'absorption' ).
%
% ARGUMENTS
% phase_shift : real 2D-array. object induced phase variations/shift
% absorption : real 2D-array. object induced absorption: log(I(x,y,z=0))
% EnergyDistancePixelsize : 3-component vector. [energy in ev,
%   sample-detector-distance in m, effective pixelsize in m]
% padding : scalar, default:0. padding factor. if 0 then no padding.
% padding_method : string, default:'symmetric'.
% verbose : boolean, default:0. print information
%
% Written by Julian Moosmann, 2018-02-05. Modified: 2018-02-06
% [int] = propagate(phase_shift, absorption, EnergyDistancePixelsize, padding, padding_method, verbose)

%% Defaults %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    phase_shift = 0.1 * phantom( 'Modified Shepp-Logan', 512 );
end
if nargin < 2
    absorption = zeros( size( phase_shift ), 'like', phase_shift );
end
if nargin < 3
    EnergyDistancePixelsize = [30e3 1 1e-6];% [eV m m]
end
if nargin < 4
    padding = 1;
end
if nargin < 5
    padding_method = 'symmetric';
end
if nargin < 6
    filt = 1;
end
if nargin < 7
    verbose = 0;
end

%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters
precision = class( phase_shift );
normalize = 1;
energy    = EnergyDistancePixelsize(1);
dist_sample_detector  = EnergyDistancePixelsize(2);
pixelsize = EnergyDistancePixelsize(3);
lambda    = E_to_lambda( energy );
frequ_prefactor = 2 * pi * lambda * dist_sample_detector / pixelsize^2;% 1/(Fresnel number) with pixelsize being the characteristic size
num_zero_crossings = frequ_prefactor / pi / 4;
im_shape = size( phase_shift );

% Padding
phase_shift = padarray( phase_shift, padding * im_shape, padding_method, 'post' );
absorption = padarray( absorption, padding * im_shape, padding_method, 'post' );
im_shape_pad = size( phase_shift );

% Fourier cooridnates
x = FrequencyVector( im_shape_pad(1), precision, normalize);
y = FrequencyVector( im_shape_pad(2), precision, normalize);
[xi, eta] = meshgrid( y, x);

% Fresnel propagation
int = abs( ifft2( filt .* exp( - 1i* frequ_prefactor * ( xi.^2 + eta.^2) / 2 ) .* fft2( exp( - absorption + 1i*phase_shift ) ) ) ).^2;

% Crop to orginial size
int = int(1:im_shape(1),1:im_shape(2));

% Print info
if verbose
    fprintf( '\nFresnel propagation:' );
    fprintf( '\n E = %g keV', energy * 1e-3 );
    fprintf( '\n lambda = %g angstrom', lambda * 1e10) 
    fprintf( '\n z = %g m', dist_sample_detector );
    fprintf( '\n pixelsize = %g micron', pixelsize * 1e6);    
    fprintf( '\n prefactor of (xi^2 + eta^2)/2 in arguments of CTFs: 2*pi*lambda*z/dx^2 = %g', frequ_prefactor );
    fprintf(' \n range of Fourier coordinates: [xi] = [%g %g], [eta] = [%g %g]', min(xi(:)), max(xi(:)), min(eta(:)), max(eta(:)) )
    fprintf( '\n shape =  [%u %u] (original), [%u %u] (padded)', im_shape, im_shape_pad );
    fprintf( '\n padding_method = %s', padding_method );      
    %fprintf('Maximum argument of the sine at [|xi| |eta|] = [1/2 1/2]: %5.3f\n',pf/4)
    fprintf( '\n CTF zero crossings - 1 (without the central one at xi = eta = 0): %5.3f', num_zero_crossings);

    fprintf( '\n' );    
    domain( phase_shift, 1, 'phase shift')
    domain( absorption, 1, 'absorption')
    domain( int, 1, 'Intensity')
    fprintf(' \n')
end
