%% Life-cell imaging, phase-contrast experiment  
% Location: at 2-BM-B, Sector 2 - Bending Magnet Beamline,APS
% Date: from 5th Feb 2012 8.00pm to 7th Feb 2012 8.00pm
%% Local Contacts:
%   Name 	XIANGHUI XIAO (micro-Tomography, Fast Tomography, Laminography)
%   Phone 	630.252.9621
%   Email 	xhxiao@aps.anl.gov
%   Name 	FRANCESCO DE CARLO (micro-Tomography)
%   Phone 	630.252.0148
%   Email 	decarlo@aps.anl.gov
% Beamline Controls and Data Acquisition: The beamline is run using Sun workstations (UNIX/Solaris). Beamline control is done though VME crates and EPICS. For diffraction experiments, data is taken using the SPEC. Some applications use the programs IDL & Igor.
% Detectors: CCD/CMOS detectors; Scintillation detectors; Energy dispersive solid state detector
%% Beamline Specs
% Source: Bending Magnet
% Monochromator Type: Multilayer monochromator
%   Energy Range:               5-30 keV
%   Resolution (ΔE/E):          1 x 10^-2
%   Flux (photons/sec):         1 x 10^12 at 17 keV
%   Beam Size (HxV), Unfocused: 25mm x 4mm
% Monochromator Type: Pink Beam
%   Energy Range:               10-30 keV
%   Flux (photons/sec):         1 x 10^14 at keV
%   Beam Size (HxV), Unfocused: 25mm x 4mm 
%% Optical Components
% Component               Distance from Source Description
% Filter assembly         23.2 m               8 filters on 2 carriers
% Hor. and vert. slits    23.5 m               25 µm reproducibility
% Vert. deflecting mirror 24.9 m               0.15 deg. plane w/ 2 coatings (Cr, Pt)
% Double multilayer mono. 27.4 m               Unfocussed
% Hor. and vert. slits    48.3 m               25 µm reproducibility 
%% Other
% 3D tomoggraphy: 1-100Hz, (2k)^2, 1500 proj, ~1min
% estimate prealignment time per sample: 2-5 min
% Scintillator: LSO 10-14 micron thicknes from ANKA, optimally 12 micron,
% LuAG at APS
% width of scintillator has to fit to depth of focus
% maximal sample translation: 5m (800mm)
% beam size after multilayer: 25mm x 4mm
% available pixelsize: 0.7, 1.4, 2.8 micron
% continous rotation
% no readout time
% prealignment time per sample: ~ 2-5min
% sample positioning motors are not encoded
% biolab available at APS
% Energy range for DMM: 5-30keV
ca
%% Experimental paramters
HutchTemp = 22.5;
sourceHor = 92e-6; %m
sourceVer = 31e-6; %m
distSampleSource = 50; %m
beamSizeHor = 25e-3; %m
beamSizeVer = 4e-3; %m
energy      = 17; %keV
wavelength  = EnergyConverter(energy);
%% Data parameter
NumProj  = 1200;
dimx     = 2048;
dimy     = 2048;
bitsize  = 16;%bit, 1byte = 8bit
NumRuns  = 7;
RunTime  = 6*3600;% seconds
TomoTime = 60; %seconds
%% Variables in BYTE!!
SizProj   = dimx*dimy*bitsize/8;
SizTomo   = NumProj*SizProj;
TotalTime = NumRuns*RunTime;
NumTomos  = NumRuns*RunTime/TomoTime;
TotalDataSiz = NumTomos*SizTomo;
% Coherence
cohFac00    = 0.61;
cohFac88    = 0.16;
cohHor00    = cohFac00*distSampleSource*wavelength/(sourceHor/2);
cohVer00    = cohFac00*distSampleSource*wavelength/(sourceVer/2);
cohHor88    = cohFac88*distSampleSource*wavelength/(sourceHor/2);
cohVer88    = cohFac88*distSampleSource*wavelength/(sourceVer/2);
%% Print
fprintf('\nDimension of projection: [%u x%u]\n',dimx,dimy)
fprintf('Data depth: %u bit (%u b)\n',bitsize,bitsize/2)
fprintf('Size of projection: %.1f Mb\n',SizProj/1024^2)
fprintf('Size of tomogram: %.1f Gb\n',SizTomo/1024^3)
fprintf('Number of runs (1 run = 6 hours = time of gastrulation): %u\n',NumRuns)
fprintf('Time for one run: %.2g h\n',RunTime/3600)
fprintf('Overall run time: %.2g h\n',TotalTime/3600)
fprintf('Number of tomograms: %g\n',NumTomos)
fprintf('Total amount of data: %.1f Gb = %.1f Tb\n',TotalDataSiz/1024^3,TotalDataSiz/1024^4)
fprintf('Coherence 08%%: horizontal: %.3g mu, vertical: %.3g mu\n',10^6*cohHor00,10^6*cohVer00)
fprintf('Coherence 88%%: horizontal: %.3g mu, vertical: %.3g mu\n',10^6*cohHor88,10^6*cohVer88)
