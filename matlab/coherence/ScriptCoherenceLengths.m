% l = R/k/a = R *lambda / 2/pi /a

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

% Coherence length in micron
l12 = (@(R,E,d) R*EnergyConverter(E)/2/pi/d*1e6);
l0 = (@(R,E,d) 3.82*R*EnergyConverter(E)/2/pi/d*1e6);

bg = (@(R,z,d) d*z/R*1e6);
