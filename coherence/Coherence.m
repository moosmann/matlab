function out = Coherence(E_keV,SourceSample_m,SourceSize_micron,SampleDetector_m,Bandwidth,PixelSize_micron)

if nargin < 1
    E_keV = 30;
end
if nargin < 2
    %SourceSample_m =  145;
    SourceSample_m = 82.7000;
end
if nargin < 3
    SourceSize_micron = sigma_to_FWHM(89.4); % PETRA III horizontal
end
if nargin < 4
    SampleDetector_m = 0.5;1;0.235;
end
if nargin < 5
   Bandwidth = 10^-4;
end
if nargin < 6
    %PixelSize_micron = 1.68 ; % micron CMOS 5 x
    PixelSize_micron = 2 * 0.64 ; % micron, CMOS 10x
end

%% Beamline parameters

% APS, 32ID, Special Operating Mode - Reduced horizontal beam size (RHB) 
% bh = 120e-6 m instead of nominal 280e-6 m
% 32ID-B: 70 m from source
% 
% http://www.aps.anl.gov/Accelerator_Systems_Division/Accelerator_Operations_Physics/SRparameters/node5.html
%
% Special Operating Mode - Reduced horizontal beam size (RHB) Lattice
% configuration: Low emittance lattice with slightly increased
% effective emittance of 3.4 nm-rad and coupling of 1%. The lattice has one
% or two locations with significantly reduced horizontal beam size: sectors
% 8 and/or 32 has horizontal beam size of 120 $\mu$m at their ID locations
% (instead of nominal 280 $\mu$m) at expence of increased divergence. All
% other radiation points are kept without changes. 

% ESRF, odd IDs: sigma = fwhm/2.355 = h: 51, v: 8.6

% DLS
% https://www.diamond.ac.uk/Science/Machine.html
%I13-2: Imaging
% source size hor 267.6 micron
% source size ver 2.9 micron
% source divergence hor 26.5 micro rad
% source divergence ver 2.8 micro rad

%I13-1: Coherence
% source size hor 307.8 micron
% source size ver 3.6  micron
% source divergence hor 18.6 micro rad
% source divergence ver 2.3 micro rad

% distance source sample 250m

%% all outpus in micron

ca

%% Constants
%h_J_s      = 6.62606896e-34;  % J * s
c      = 299792458; % m / s
%q      = 1.60217733e-16;

h_eV_s = 4.135667516e-15; % eV * s

%% MAIN

lambda = EnergyConverter(E_keV);
R = SourceSample_m;
s = SourceSize_micron * 1e-6;
z = SampleDetector_m;
PixelSize = PixelSize_micron * 1e-6; % m

%% Coherence length in micron
% circular source, Airy function drop to .88 or 0
lc88 = (@(lambda,R,s)       R * lambda / (2 * pi * s) * 1e6 );
lc00 = (@(lambda,R,s) 3.82 * R * lambda / (2 * pi * s) * 1e6 );
% Gaussian source, drops to 0.5
lg50 = (@(lambda,R,s)  4 * log(2) * lc88(lambda,R,s) );
lg50approx = (@(lambda,R,s)  R * lambda / (2 * s) * 1e6 );

ll = @ ( Bandwidth, E_keV ) h_eV_s * c / ( Bandwidth * E_keV * 1e3 ) * 1e6;  % 2 pi v / domega
llt = @ (R, lcohl) sqrt( 2 * R * lcohl/1e6 - (lcohl/1e6) ^ 2) * 1e6;
% source blurring
bg = (@(R,z,s) z * s / R * 1e6);

%% Out struct
lcohl = ll(Bandwidth,E_keV);
% coherenece
out.coherence.circular88 = lc88(lambda,R,s);
out.coherence.circular00 = lc00(lambda,R,s);
out.coherence.gaussian50 = lg50(lambda,R,s);
out.coherence.gaussian50approx = lg50approx(lambda,R,s);
out.coherence.longitudinal = lcohl;
out.coherence.longToTrans = llt(R,lcohl);
% frequency, angular
out.omega.zero = E_keV / h_eV_s * 2 * pi;
out.omega.delta = Bandwidth * E_keV / h_eV_s * 2 * pi;
% blur
out.blur.geometric = bg(R,z,s);
% spatial frequency cut-off
out.frequency.blurOverPixel = PixelSize_micron / out.blur.geometric;
% distance of points that interfere due to Guigay's autocorrelation
out.distance.interferenceByAutocorrelation = lambda * z / 2 / (PixelSize) * 1e6;

%% Print %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\nCohrence lengths in microns: \n\n')
disp(out.coherence)

fprintf('\n angular frequency: %g +/- %g s', out.omega.zero, out.omega.delta)

fprintf( '\n energy: %g keV', E_keV )
fprintf( '\n propagation distance: %g m', SampleDetector_m )
fprintf( '\n effective pixel size: %g micron', PixelSize_micron )

fprintf('\n geometric blur: %g micron', out.blur.geometric)
fprintf('\n geom. blur / pixelsize: %g micron', out.blur.geometric / PixelSize_micron)

fprintf('\n frequency cut-off by geometric blur (pixel size / geomtric blur):\n  %g', out.frequency.blurOverPixel)

fprintf('\n distance between sample points that interfere due to autocorrelation in micron:\n  %g micron', out.distance.interferenceByAutocorrelation)

%% Plot
zRange = 0:0.005:2;
b = bg(R,zRange,s);
[~, zmaxind] = min( abs( b - PixelSize_micron ) );
zmax = zRange(zmaxind);


figure( 'Name', 'geometric blur vs propagation distance')
plot( zRange, [b; repmat( PixelSize_micron, numel( zRange))] )
xlabel( 'propagation distance / mi' )
ylabel( 'spatial resolution / micron' )
legend( 'geometric blur', 'effective pixelsize' )
axis tight

fprintf( ' \n distance where blur equals eff pixel size: %g m', zmax )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\n' )