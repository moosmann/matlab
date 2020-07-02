%% BM5 Characteristics
% Source :
% Bending magnet
% Magnetic field: 0.82 Tesla
% Radius of curvature: 24.593 m
% Critical energy: 19.87 keV
% X-ray source characteristics (at 8 mrad)
% H. (FWHM) size: 270 μm
% V. (FWHM) size: 80 μm
% H. (FWHM) divergence: 2.4 mrad
% V. (FWHM) divergence: 180 μrad
% Power:120 W/mrad, 1.35 W/mm2
% Max. flux: 2.7*10^13 ph./s/mrad2/0.1 %BW
% Optics :
% Double crystal Si(111) monochromator
% Energy range 6-60 keV
% Energy resolution : E/E~ 2.1 10-4
% Distance from source : 27.22 m
% Flux at 25 keV and 200 mA:
% • With flat crystals: 1.6*10^10 ph/s in 1 mm x 1 mm
% • Using sagittal crystal: 2.2*10^12 ph/s in 1 mm x 0.3 mm
% Double-multilayer monochromator (6 keV- 30 keV)
% Energy range 6-30 keV
% Energy resolution: dE/E~ 3.7 10-2
% Distance from source: 28.4 m
%% APS
% Beamline Specs
% Source: Bending Magnet
% Monochromator Type: Multilayer monochromator
%   Energy Range: 5-30 keV
%   Resolution (ΔE/E): 1*10^-2
%   Flux (photons/sec): 1*10^12 @17 keV
%   Beam Size (HxV), Unfocused: 25mm x 4mm
% Monochromator Type: Pink Beam
%   Energy Range: 10-30 keV
%   Flux (photons/sec): 1*10^14 @ keV
%   Beam Size (HxV), Unfocused: 25mm x 4mm 

% BM05: Flux at ~25 keV, 200 mA (bending magnet), 30m, and 10^-2 bandwidth (multilayer):
flux = 10*140*10^9; % photons/s/(1mm)^2/0.1%
% total exposure time for one tomogram
tet = 120; % s
% area illuminated by the beam
ia = 4; % (mm)^2
% Number of photons absorbed in the container
noap = 0.5*flux*tet*ia;
% Energy deposited due to absorption in the container
ed = noap*25e3; % eV
disp(flux)
disp(noap)
disp(ed)
