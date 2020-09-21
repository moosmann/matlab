function out = Coherence(E_keV,SourceSample_m,SourceSize_micron,SampleDetector_m,Bandwidth,PixelSize_micron)

if nargin < 1
    E_keV = 30;
end
if nargin < 2
    %SourceSample_m =  145;
    SourceSample_m = 87;
end
if nargin < 3
    SourceSize_micron = sigma_to_FWHM( 90 );% P05 low beta section
    %sigma_to_FWHM(89.4); % PETRA III horizontal
end
if nargin < 4
    SampleDetector_m = 1.4;
end
if nargin < 5
   Bandwidth = 10^-4;
end
if nargin < 6
    %PixelSize_micron = 1.68 ; % micron CMOS 5 x
    PixelSize_micron = 2 ; % micron, CMOS 10x
end

%% Beamline parameters

% P05
% The Imaging Beamline P05 (IBL) is dedicated to full-field imaging
% techniques. IBL is located on the low-β sector 4 of the PETRA III hall
% and shares this sector with the Hard X-ray Micro/Nano-Probe Beamline P06.
% The insertion device is a 2 m long U29 undulator of a 5 mrad canted
% undulator pair. The designed size (rms) of the source at 10 keV is 36 μm
% x 6.1 μm with a divergence (rms) of 28 x 4.0 μrad2. Apart from the
% undulator, the beamline frontend contains a power slit system and filters
% (4 mm glassy carbon and 300 μm diamond window with 50 μm Cu). The
% beamline is designed to operate in an energy range from 5 - 50 keV. 
%
% IBL
% The microtomography endstation (EH2) provides imaging techniques with
% spatial resolutions down to 0.9 μm. The sample is located ~87 m from the
% source, at which point the beam has diverged to approximately 7 x 2.5 mm2
% (DCM, channel cut) or 11 x 11 mm2 (DMM) (FWHM).   


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
font_size = 18;
line_width = 6;

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
fprintf('\nCoherence lengths in microns: \n\n')
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
zRange = 0:0.005:1.4;
b = bg(R,zRange,s);
[~, zmaxind] = min( abs( b - PixelSize_micron ) );
zmax = zRange(zmaxind);

fig_path = '/asap3/petra3/gpfs/common/p05/jm/pictures/';
h1 = figure( 'Name', 'Coherence: geometric blur vs propagation distance');
p = plot( zRange, [b; repmat( PixelSize_micron, numel( zRange))] );
xlabel( 'propagation distance / m' )
ylabel( 'spatial resolution / micron' )
legend( 'geometric blur', 'effective pixelsize', 'Location', 'northwest' )
title( sprintf( 'energy: %g keV', E_keV) )
set( p ,'LineWidth', line_width )
ax = gca;
ax.FontSize = font_size; 
axis tight
saveas( h1, sprintf( '%s%s.png', fig_path, regexprep( h1.Name, '\ |:', '_') ) );

fprintf( ' \n distance where blur equals eff pixel size: %g m', zmax )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\n' )