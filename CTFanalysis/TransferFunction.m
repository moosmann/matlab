function [transferFunction] = TransferFunction(Method,imSize,EnergyDistancePixelsize,RegPar,BinaryFilterThreshold,outputPrecision)
% Transfer function related to 'PhaseFilter'
% 
% Written by Julian Moosmann, last modification: 2014-10-30
%

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    Method = 'tie';
end
if nargin < 2
    imSize = [1024 1024];
end
if nargin < 3
    EnergyDistancePixelsize = [20 0.945 .75e-6];
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

%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call 'PhaseFilterDual' if duality is assumed. Then, 'regPar' is
% interpreted as -log10 of duality factor epsilon.

%% Standard phase retrieval %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameter
Energy    = EnergyDistancePixelsize(1);
Distance  = EnergyDistancePixelsize(2);
Pixelsize = EnergyDistancePixelsize(3);
% wave length
lambda    = 6.62606896e-34*299792458/(Energy*1.60217733e-16);
ArgPrefac = 2*pi*lambda*Distance/Pixelsize^2;

%% Fourier coordinates
% 1D
xi  = FrequencyVector(imSize(2),outputPrecision,1);
eta = FrequencyVector(imSize(1),outputPrecision,1);
% 2D
[sinArg, sinxiquad]   = meshgrid(xi,eta);
% Function on 2D
sinArg = ArgPrefac*(sinArg.^2 + sinxiquad.^2)/2;

%% Filter 
switch lower(Method)
    case 'tie'
        transferFunction  = 2*(sinArg + 10^-RegPar);
    case 'ctf'
        sinxiquad   = sin(sinArg);
        transferFunction  = 2*sign(sinxiquad).*(abs(sinxiquad)+10^-RegPar);        
    case {'ctfhalfsine','ctffirsthalfsine','halfsine','firsthalfsine'}
        sinxiquad   = sin(sinArg);
        transferFunction  = 2*sign(sinxiquad).*(abs(sinxiquad)+10^-RegPar);
        transferFunction ( sinArg >= pi ) = 0;
    case {'qp','pctf','quasi','quasiparticle','quasiparticles'}
        sinxiquad   = sin(sinArg);
        transferFunction  = 2*sign(sinxiquad).*(abs(sinxiquad)+10^-RegPar);
        transferFunction ( sinArg > pi/2  &  abs(sinxiquad) < BinaryFilterThreshold) = 0;
    case {'qphalfsine','pctfhalfsine','pctffirsthalfsine','quasihalfsine','quasifirsthalfsine'}
        sinxiquad   = sin(sinArg);
        transferFunction  = 2*sign(sinxiquad).*(abs(sinxiquad)+10^-RegPar);
        transferFunction ( sinArg > pi/2  &  abs(sinxiquad) < BinaryFilterThreshold) = 0;
        transferFunction ( sinArg >= pi ) = 0;
    case {'qp2','quasi2','pctf2'}
        sinxiquad   = sin(sinArg);
        transferFunction  = 2*sign(sinxiquad).*(abs(sinxiquad)+10^-RegPar);
        mask = sinArg > pi/2  &  abs(sinxiquad) < BinaryFilterThreshold;
        transferFunction ( mask ) = bsxfun(@(a,b) a(b),sign(transferFunction )/(2*(BinaryFilterThreshold+10^-RegPar)),mask);
end

