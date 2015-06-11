function [phaseFilter,phaseAppendix] = PhaseFilter(Method,imSize,EnergyDistancePixelsize,RegPar,BinaryFilterThreshold,outputPrecision)
% Fourier space filter for single-distance phase retrieval. 
%
% Given a normalized intensity map I_z (i.e. flat- and dark-field
% corrected), the intensity contrast is defined as g_z = I_z - 1. 
% Phase retrieval: phi = real( ifft2( PhaseFilter() .* fft2( g_z ) ) ). No
% fftshift is needed. Optional output is a string consisting of the
% parameters used. 
%
% Arguments:
%
% Method : string. Default: 'tie'.
%   Variants:
%   'tie': Linearized transport of intensity equation. Essentially the
%   inversion of the Laplacian, sometimes also referred to as Paganing
%   or Bronnikov phase retrieval.
%   'ctf': Inversion of the contrast transfer function in the pure phase
%   case. Is essentially the inversion of a sine.
%   'ctfhalfsine': Same as 'ctf' up to the first half period of the sine.
%   'qp': Quasiparticle version of 'ctf', i.e. cropping of frequency bands
%   centered around the positions of the zero crossing of the Fourier
%   transform of g_z. See papers: http://dx.doi.org/10.1364/OE.19.025881
%   and http://dx.doi.org/10.1364/OE.19.012066.
%   'qphalfsine': use only first half period of the sine of quasiparticle
%   version of 'ctf'.
%
% imSize : 1x2-vector. Default [1024 1024]. Size of output filter.
%
% EnergyDistancePixelsize : 1x3-vector. Default: [20 0.945 .75e-6]. Energy
% in keV, Distance in m, Pixelsize in m.
%
% RegPar : scalar. Default: 2.5. Phase retrieval is regularised according:
% 1/func(x) -> 1/(func(x)+10^(-RegPar)), for details of the placeholder
% function func see below. Thus, the regularization parameter RegPar is the
% negative of the decadic logartihm of the constant which is added to the
% denominator in order to regularize the singularity at zero frequency.
% Typical values are between 1.5 and 3.5 depending on energy, residual
% absorption, etc.
%
% BinaryFilterThreshold : scalar. Default: 0.1. Parameter for Quasiparticle
% phase retrieval defining the width of the rings which are cropped around
% the zero crossings of the CTF denominator (in Fourier space). Typical
% values are between 0.01 and 0.1, where large values yields results
% similiar to 'tie' phase retrieval.
%
% outputPrecision : string. 'single' (default), or 'double'.
% 
% Written by Julian Moosmann, last modification: 2015-05-27
%
% [phaseFilter,phaseAppendix] = PhaseFilter(Method,imSize,EnergyDistancePixelsize,RegPar,BinaryFilterThreshold,outputPrecision)

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
strInd = strfind(lower(Method),'dual');
switch ~isempty(strInd)
%% Standard phase retrieval %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 0
%% Parameter
Energy    = EnergyDistancePixelsize(1);
Distance  = EnergyDistancePixelsize(2);
Pixelsize = EnergyDistancePixelsize(3);
% wave length
lambda    = 6.62606896e-34*299792458/(Energy*1.60217733e-16);
ArgPrefac = pi*lambda*Distance/Pixelsize^2;

%% Fourier coordinates
% 1D
xi  = FrequencyVector(imSize(2),outputPrecision,1);
eta = FrequencyVector(imSize(1),outputPrecision,1);
% 2D
[sinArg, sinxiquad]   = meshgrid(xi,eta);
% Function on 2D
sinArg = ArgPrefac*(sinArg.^2 + sinxiquad.^2);

%% Filter 
switch lower(Method)
    case 'tie'
        phaseFilter = 1/2./(sinArg + 10^-RegPar);
        phaseAppendix = sprintf('tie_regPar%3.2f',RegPar);
    case 'ctf'
        sinxiquad   = sin(sinArg);
        phaseFilter = 1/2*sign(sinxiquad)./(abs(sinxiquad)+10^-RegPar);        
        phaseAppendix = sprintf('ctf_regPar%3.2f',RegPar);
    case {'ctfhalfsine','ctffirsthalfsine','halfsine','firsthalfsine'}
        sinxiquad   = sin(sinArg);
        phaseFilter = 1/2*sign(sinxiquad)./(abs(sinxiquad)+10^-RegPar);
        phaseFilter( sinArg >= pi ) = 0;
        phaseAppendix = sprintf('ctfHalfSine_regPar%3.2f',RegPar);
    case {'qp','pctf','quasi','quasiparticle','quasiparticles'}
        sinxiquad   = sin(sinArg);
        phaseFilter = 1/2*sign(sinxiquad)./(abs(sinxiquad)+10^-RegPar);
        phaseFilter( sinArg > pi/2  &  abs(sinxiquad) < BinaryFilterThreshold) = 0;
        phaseAppendix = sprintf('qp_regPar%3.2f_binFilt%3.3f',RegPar,BinaryFilterThreshold);
    case {'qphalfsine','pctfhalfsine','pctffirsthalfsine','quasihalfsine','quasifirsthalfsine'}
        sinxiquad   = sin(sinArg);
        phaseFilter = 1/2*sign(sinxiquad)./(abs(sinxiquad)+10^-RegPar);
        phaseFilter( sinArg > pi/2  &  abs(sinxiquad) < BinaryFilterThreshold) = 0;
        phaseFilter( sinArg >= pi ) = 0;
        phaseAppendix = sprintf('qpHalfSine_regPar%3.2f_binFilt%3.3f',RegPar,BinaryFilterThreshold);
    case {'qp2','quasi2','pctf2'}
        sinxiquad   = sin(sinArg);
        phaseFilter = 1/2*sign(sinxiquad)./(abs(sinxiquad)+10^-RegPar);
        mask = sinArg > pi/2  &  abs(sinxiquad) < BinaryFilterThreshold;
        phaseFilter( mask ) = bsxfun(@(a,b) a(b),sign(phaseFilter)/(2*(BinaryFilterThreshold+10^-RegPar)),mask);
        phaseAppendix = sprintf('qp2_regPar%3.2f_binFilt%3.3f',RegPar,BinaryFilterThreshold);
end

% Restore zero frequency component
phaseFilter(1) = 1/2*10^RegPar;

% Replace dots by p
phaseAppendix = regexprep(phaseAppendix,'\.','p');

%% Dual phase retrieval %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 1
        % remove string 'dual' from Variable 'Method' to use it as input
        % for 'PhaseFilterDual'
        Method(strInd+(0:3)) = [];
        % Call 'PhaseFilterDual'
        [phaseFilter,phaseAppendix] = PhaseFilterDual(Method,imSize,EnergyDistancePixelsize,10^-RegPar,BinaryFilterThreshold,outputPrecision);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function k = FrequencyVector(SizeOfVector,Precision,Normalize)
% Frequency vector for MATLAB's implementation of Fourier space
% discretization (without 'fftshift') ranging from 0 to
% floor(('SizeOfVector'-1)/2) and -floor('SizeOfVector'/2) to -1.
%
% SizeOfVector: integer.
%
% Precsision: string, default: 'single'. Available: 'single', or 'double'.
%
% Normalize: boolean, default: 1. Output frequencies normalized by size of
% vector.
%
% Written by Julian Moosmann, last version 2013-11-12.
%
%k = FrequencyVector(SizeOfVector,Precision,Normalize)

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    Precision = 'single';
end
if nargin < 3
    Normalize = 1;
end

%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Precision
SizeOfVector = eval(sprintf('%s(%u)',Precision,SizeOfVector));

% Create vector
k = [0:1:ceil( (SizeOfVector-1)/2 )  -floor( (SizeOfVector-1)/2):1:-1] /(1 -Normalize*(1 - 1*SizeOfVector ) );