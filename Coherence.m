function out = Coherence(E_keV,SourceSample_m,SourceSize_micron,SampleDetector_m,Bandwidth,PixelSize_micron)

if nargin < 1
    E_keV = 30;
end
if nargin < 2
    SourceSample_m =  145;
end
if nargin < 3
    SourceSize_micron = SigmaToFWHM(90);
end
if nargin < 4
    SampleDetector_m = 2;
end
if nargin < 5;
   Bandwidth = 10^-4;
end
if nargin < 6
    PixelSize_micron = 1.6 ; % micron
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


%% all outpus in micron


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

fprintf('Angular frequency in seconds: \n\n')
disp(out.omega)

fprintf('Geometric blur in micron: \n\n')
disp(out.blur)

fprintf('Blur-induced over resolution-induced frequency cut-off: \n\n')
disp(out.frequency)

fprintf('Distance between sample points that interfere due to autocorrelation in micron: \n\n')
disp(out.distance)

