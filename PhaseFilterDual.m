function [phaseFilter,phaseAppendix] = PhaseFilterDual(Method,imSize,EnergyDistancePixelsize,Epsilon,BinaryFilterThreshold,outputPrecision)
% Phase retrieval filter in Fourier space: 2D filter to be multiplied by
% the Fourier transformed intensity to retrieve the Fourier transform of
% the phase.
%
% For details see help of 'PhaseFilter'. In contrast to 'PhaseFilter' this
% function assume a duality (proportionality) between attenuation B and
% phase to hold: B = -Epsilon * phase
%
% Written by Julian Moosmann, last version: 2014-01-15
%
% [phaseFilter,phaseAppendix] = PhaseFilterDual(Method,imSize,EnergyDistancePixelsize,Epsilon,BinaryFilterThreshold,outputPrecision) 

% Refractive index of water H2O at 20 and 30 keV
% Energy/eV  Delta,          Beta            Beta/Delta  -log10(Beta/Delta)
% 20000      5.76620721E-07  3.46201373E-10  6.0040e-04  3.2216
% 30000      2.56114134E-07  1.06431752E-10  4.1556e-04  3.3814

%% Default arguments
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
    Epsilon = 2.5;
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
lambda    = 6.62606896e-34*299792458/(Energy*1.60217733e-16);
ArgPrefac = 2*pi*lambda*Distance/Pixelsize^2;
%% Fourier coordinates.
% 1D
%xi  = -1/2:1/imSize(2):1/2-1/imSize(2);
xi  = FrequencyVector(imSize(2),outputPrecision,1);
%eta = -1/2:1/imSize(1):1/2-1/imSize(1);
eta = FrequencyVector(imSize(1),outputPrecision,1);
% 2D
[sinArg, sinxiquad]   = meshgrid(xi,eta);
% Function on 2D
sinArg = ArgPrefac*(sinArg.^2 + sinxiquad.^2)/2 + atan(Epsilon);
%% Filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch lower(Method)
    case 'tie'
        phaseFilter = 1/sqrt(1+Epsilon^2)/2./sinArg;
        phaseAppendix = sprintf('tieDual_mlogeps%3.2f',-log10(Epsilon));
    case 'ctfold'
        sinxiquad   = sin(sinArg);
        phaseFilter = 1/sqrt(1+Epsilon^2)/2*sign(sinxiquad)./abs(sinxiquad);        
        phaseAppendix = sprintf('ctfDual_mlogeps%3.2f',-log10(Epsilon));
    case 'ctf'
        sinxiquad   = sin(sinArg);
        phaseFilter = 1/sqrt(1+Epsilon^2)/2*sign(sinxiquad)./(abs(sinxiquad) + Epsilon);        
        phaseAppendix = sprintf('ctfDualReg_mlogeps%3.2f',-log10(Epsilon));       
    case {'ctfhalfsine','ctffirsthalfsine','halfsine','firsthalfsine'}
        sinxiquad   = sin(sinArg);
        phaseFilter = 1/sqrt(1+Epsilon^2)/2*sign(sinxiquad)./abs(sinxiquad);
        phaseFilter( sinArg >= pi ) = 0;
        phaseAppendix = sprintf('ctfHalfSineDual_mlogeps%3.2f',-log10(Epsilon));
    case {'qphalfsine','pctfhalfsine','pctffirsthalfsine','quasihalfsine','quasifirsthalfsine'}
        sinxiquad   = sin(sinArg);
        phaseFilter = 1/sqrt(1+Epsilon^2)/2*sign(sinxiquad)./abs(sinxiquad);
        phaseFilter( sinArg > pi/2  &  abs(sinxiquad) < BinaryFilterThreshold) = 0;
        phaseFilter( sinArg >= pi ) = 0;
        phaseAppendix = sprintf('quasiHalfSineDual_mlogeps%3.2f_binFilt%3.3f',-log10(Epsilon),BinaryFilterThreshold);    
    case {'qp','pctf','quasi','quasiparticle','quasiparticles'}        
        sinxiquad   = sin(sinArg);
        phaseFilter = 1/sqrt(1+Epsilon^2)/2*sign(sinxiquad)./abs(sinxiquad);
        phaseFilter( sinArg > pi/2 &  abs(sinxiquad) < BinaryFilterThreshold) = 0;   
        phaseAppendix = sprintf('qpDual_mlogeps%3.2f_binFilt%3.3f',-log10(Epsilon),BinaryFilterThreshold);    
    case {'qp2','quasi2','pctf2'}
        sinxiquad   = sin(sinArg);
        phaseFilter = 1/sqrt(1+Epsilon^2)/2*sign(sinxiquad)./abs(sinxiquad);
        mask = sinArg > pi/2  &  abs(sinxiquad) < BinaryFilterThreshold;
        phaseFilter( mask ) = bsxfun(@(a,b) a(b),sign(phaseFilter)/(2*(BinaryFilterThreshold)),mask);
        phaseAppendix = sprintf('qp2Dual_mlogeps%3.2f_binFilt%3.3f',-log10(Epsilon),BinaryFilterThreshold);
end

% Restore zero frequency component
phaseFilter(1) = 1/(2*Epsilon);

% Replace dots by p
phaseAppendix = regexprep(phaseAppendix,'\.','p');
