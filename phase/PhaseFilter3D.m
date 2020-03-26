function [phaseFilter,phaseAppendix] = PhaseFilter3D(Method,vol_shape,EnergyDistancePixelsize,RegPar,BinaryFilterThreshold,outputPrecision)
% Phase filter for reconstruction of the real refractive index decrement
% 'delta' on a tomographic volume. Parameters as in 'PhaseFilter'.
%
% Written by Julian Moosmann, last version: 2013-11-13
%
% [phaseFilter,phaseAppendix] = PhaseFilter3D(Method,vol_shape,EnergyDistancePixelsize,RegPar,BinaryFilterThreshold,outputPrecision)

if nargin < 1
    Method = 'tie';
end
if nargin < 2
    vol_shape = [1024 1024 1024];
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

% %% Fourier coordinates.
% % 1D
% xi  = FrequencyVector(vol_shape(2),outputPrecision,1);
% eta = FrequencyVector(vol_shape(1),outputPrecision,1);
% zeta = FrequencyVector(vol_shape(3),outputPrecision,1);
% % 3D
% [sinArg, sinxiquad, zeta]  = meshgrid(xi,eta,zeta);
% % Function on 2D
% sinArg = ArgPrefac*(sinArg.^2 + sinxiquad.^2 + zeta.^2)/2;

%% Fourier space variables
% 1D
xi1(:,1,1) = FrequencyVector(vol_shape(1),outputPrecision,1);
xi2(1,:,1) = FrequencyVector(vol_shape(2),outputPrecision,1);
xi3(1,1,:) = FrequencyVector(vol_shape(3),outputPrecision,1);
% 3D
sinArg = ArgPrefac*(repmat(xi1,[1 vol_shape(2) vol_shape(3)]).^2 + repmat(xi2,[vol_shape(1) 1 vol_shape(3)]).^2 + repmat(xi3,[vol_shape(1) vol_shape(2) 1]).^2)/2;

%% Filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch lower(Method)
    case 'tie'
        phaseFilter = 1/2./(sinArg + 10^-RegPar);
        phaseAppendix = sprintf('tie_regPar%3.2f',RegPar);
    case 'ctf'
        sinArg   = sin(sinArg);
        phaseFilter = 1/2*sign(sinArg)./(abs(sinArg)+10^-RegPar);
        phaseAppendix = sprintf('ctf_regPar%3.2f',RegPar);
    case {'ctfhalfsine','ctffirsthalfsine','halfsine','firsthalfsine'}
        phaseFilter   = sin(sinArg);
        phaseFilter = 1/2*sign(phaseFilter)./(abs(phaseFilter)+10^-RegPar);
        phaseFilter( sinArg >= pi ) = 0;
        phaseAppendix = sprintf('ctfHalfSine_regPar%3.2f',RegPar);
    case {'qphalfsine','quasihalfsine','quasifirsthalfsine'}
        sinxiquad   = sin(sinArg);
        phaseFilter = 1/2*sign(sinxiquad)./(abs(sinxiquad)+10^-RegPar);
        phaseFilter( sinArg > pi/2  &  abs(sinxiquad) < BinaryFilterThreshold) = 0;
        phaseFilter( sinArg >= pi ) = 0;
        phaseAppendix = sprintf('qpHalfSine_regPar%3.2f_binFilt%3.3f',RegPar,BinaryFilterThreshold);
    case {'qp','pctf','quasi','quasiparticle','quasiparticles'}
        sinxiquad   = sin(sinArg);
        phaseFilter = 1/2*sign(sinxiquad)./(abs(sinxiquad)+10^-RegPar);
        phaseFilter( sinArg > pi/2  &  abs(sinxiquad) < BinaryFilterThreshold) = 0;
        %         phaseFilter = 1/2*sign(sin(sinArg))./(abs(sin(sinArg))+10^-RegPar);
        %         phaseFilter( sinArg > pi/2  &  abs(sin(sinArg)) < BinaryFilterThreshold) = 0;
        phaseAppendix = sprintf('qp_regPar%3.2f_binFilt%3.3f',RegPar,BinaryFilterThreshold);
    case {'qp2','qpnew','quasinew','quasi2','pctf2','sinesupression'}
        sinxiquad   = sin(sinArg);
        phaseFilter = 1/2*sign(sinxiquad)./(abs(sinxiquad)+10^-RegPar);
        mask = sinArg > pi/2  &  abs(sinxiquad) < BinaryFilterThreshold;
        phaseFilter( mask ) = bsxfun(@(a,b) a(b),sign(phaseFilter)*1/(2*(BinaryFilterThreshold+10^-RegPar)),mask);
        phaseAppendix = sprintf('qp2_regPar%3.2f_binFilt%3.3f',RegPar,BinaryFilterThreshold);
end

%% Set zero frequency back
phaseFilter(1) = 1/2*10^RegPar;

%% Replace dots by p
phaseAppendix = regexprep(phaseAppendix,'\.','p');
