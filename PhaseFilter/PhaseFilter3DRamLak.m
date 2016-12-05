function [phaseFilter,phaseAppendix] = PhaseFilter3DRamLak(Method,imSize,EnergyDistancePixelsize,RegPar,BinaryFilterThreshold,outputPrecision)
% Phase filter for reconstruction of the real refractive index decrement
% 'delta' on a tomographic volume. Parameters as in 'PhaseFilter'.
%
% Written by Julian Moosmann, last version: 2016-12-05
%
% [phaseFilter,phaseAppendix] = PhaseFilter3D(Method,imSize,EnergyDistancePixelsize,RegPar,BinaryFilterThreshold,outputPrecision)

% !!!! PHASE RETRIEVAL REGULARIZATION FOR CTF MODELS IS BUGGY. !!!!!!!!!!!

if nargin < 1
    Method = 'tie';
end
if nargin < 2
    imSize = [1024 1024 1024];
end
if nargin < 3
    EnergyDistancePixelsize = [30e3 .4 1e-6];
end
if nargin < 4
    RegPar = 2.5;
end
if nargin < 5
    BinaryFilterThreshold = 0.1;
end
if nargin < 6
    outputPrecision = 'single';
end

%% Parameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Energy    = EnergyDistancePixelsize(1);
Distance  = EnergyDistancePixelsize(2);
Pixelsize = EnergyDistancePixelsize(3);
lambda    = 6.62606896e-34*299792458/(Energy/1000*1.60217733e-16);
% Prefactor needed for TIE and CTF retrieval.
ArgPrefac = 2*pi*lambda*Distance/Pixelsize^2;

%% Fourier coordinates.
% 1D
%xi  = -1/2:1/imSize(2):1/2-1/imSize(2);
 xi  = FrequencyVector(imSize(2),outputPrecision,1);
%eta = -1/2:1/imSize(1):1/2-1/imSize(1);
eta = FrequencyVector(imSize(1),outputPrecision,1);
%zeta =-1/2:1/imSize(3):1/2-1/imSize(3);
zeta = FrequencyVector(imSize(3),outputPrecision,1);
% 2D
[xi, eta, zeta]  = meshgrid(xi,eta,zeta);
% Function on 2D
sinArg = ArgPrefac*(xi.^2 + eta.^2 + zeta.^2)/2;

%% Filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch lower(Method)
    case 'tie'
        phaseFilter = 1/2./(sinArg + 10^-RegPar);
        phaseAppendix = sprintf('tie_regPar%3.2f',RegPar);
    case 'ctf'
        sinxiquad   = sin(sinArg);
        phaseFilter = 1./(2*sign(sinxiquad).*(abs(sinxiquad))+10^-RegPar);
        phaseAppendix = sprintf('ctf_regPar%3.2f',RegPar);
    case {'ctfhalfsine','ctffirsthalfsine','halfsine','firsthalfsine'}
        sinxiquad   = sin(sinArg);
        phaseFilter = 1./(2*sign(sinxiquad).*(abs(sinxiquad))+10^-RegPar);
        phaseFilter( sinArg >= pi ) = 0;
        phaseAppendix = sprintf('ctfHalfSine_regPar%3.2f',RegPar);
    case {'pctf','quasi','quasiparticle','quasiparticles'}
        sinxiquad   = sin(sinArg);
        phaseFilter = 1./(2*sign(sinxiquad).*(abs(sinxiquad))+10^-RegPar);
        phaseFilter( sinArg > pi/2  &  abs(sinxiquad) < BinaryFilterThreshold) = 0;
        phaseAppendix = sprintf('quasi_regPar%3.2f_binFilt%3.3f',RegPar,BinaryFilterThreshold);
    case {'pctfhalfsine','pctffirsthalfsine','quasihalfsine','quasifirsthalfsine'}
        sinxiquad   = sin(sinArg);
        phaseFilter = 1./(2*sign(sinxiquad).*(abs(sinxiquad))+10^-RegPar);
        phaseFilter( sinArg > pi/2  &  abs(sinxiquad) < BinaryFilterThreshold) = 0;
        phaseFilter( sinArg >= pi ) = 0;
        phaseAppendix = sprintf('quasiHalfSine_regPar%3.2f_binFilt%3.3f',RegPar,BinaryFilterThreshold);
    case {'quasinew','quasiquasi','quasi2','pctf2','sinesupression'}
        sinxiquad   = sin(sinArg);
        phaseFilter = 1./(2*sign(sinxiquad).*(abs(sinxiquad))+10^-RegPar);
        mask = sinArg > pi/2  &  abs(sinxiquad) < BinaryFilterThreshold;
        phaseFilter( mask ) = bsxfun(@(a,b) a(b),sign(phaseFilter)*1/(2*(BinaryFilterThreshold+0*10^-RegPar)),mask);
        phaseAppendix = sprintf('quasiNew_regPar%3.2f_binFilt%3.3f',RegPar,BinaryFilterThreshold);
end

%% Ram-Lak
%RamLak = @(VolSize) repmat(sqrt(repmat(FrequencyVector(VolSize(2)),[VolSize(1) 1]).^2 + repmat(FrequencyVector(VolSize(1))',[1 VolSize(2)]).^2),[1 1 VolSize(3)]);
phaseFilter = phaseFilter.* sqrt(xi.^2 + eta.^2);
%% String
% Replace dots by p
phaseAppendix = regexprep(['ramLak_' phaseAppendix],'\.','p');
