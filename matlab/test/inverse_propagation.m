
N = 1000;

% phase object
phase = zeros(N,1);
phase(ceil(N/4):N-floor(N/4)) = 1;
phase = FilterBlur( phase );

% Parameters.
energy      = 20; % keV
lambda      = EnergyConverter(energy); % in m
distance    = 0.5; % m
pixelsize   = 1e-6; % in m

% Fourier cooridnates.
xi = fftshift(-1/2:1/N:1/2-1/N)';

% Propagator.
fprop = exp( -1i * pi * lambda * distance / pixelsize^2 * xi.^2);

% inverse propagator
fprop_inv = 1 ./ fprop;

% transmission
trans = exp( 1i * phase );

% FT of transmission
ftrans = fft( trans, [], 1 );

% wave field
w = ifft ( fprop .* ftrans, [], 1 );

% intensity of forward propagated wave field
int = abs( w ) .^2;

% backpropagated wave field
wb = ifft( fprop_inv .* fft( w ) );

% intensity of backpropagated wave field
intb = abs( wb ) .^2;

% backpropagated phase
phaseb = imag( log( wb ) );

plot( phase - phaseb )
