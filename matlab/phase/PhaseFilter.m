function [fourier_filter, parameter_string] = PhaseFilter(method, filter_size, energy_distance_pixelsize, regularization_parameter, binary_filter_threshold, frequency_cutoff, precision)
% Fourier space filter for Fourier-transform-based algebraic
% single-distance phase retrieval.
%
% Phase retrieval:
%
%   pha = real( ifft2( fourier_filter(ARG) .* fft2( g_z ) ) ),
% where g_z is the intensity contrast defined as g_z = I_z - 1 with the
% normalized intensity map I_z (i.e. flat- and dark-field corrected).
%
% OUTPUT:
%
% fourier_filter: 2D array, Fourier space filter for algebraic
%   fourier-based phase retrieval. % no fftshift required.
% parameter_sring: string, optional, unique string summarizing the
%   parameters used.
% 
% ARGUMENTS:
%
% method : string. Default: 'tie'.
%   Variants:
%   'tie' : Linearized transport of intensity equation. Amounts to the
%     inversion of the Laplacian, also referred to as Paganin
%     or Bronnikov phase retrieval.
%   'ict' : Intensity contrast transfer assuming weak phase and weak
%   attenaution variations, but not a weak absolute attenuation. 
%       https://doi.org/10.1364/OL.530330
%   'ctf' : Inversion of the contrast transfer function for the pure phase
%     case. Amounts to the inversion of the CTF sinusoidal prefactor.
%   'ctfhalfsine': Same as 'ctf' up to the first half period of the sine.
%   'qp' : Quasiparticle version of 'ctf', i.e. cropping of frequency bands
%     centered around the positions of the zero crossing of the Fourier
%     transform of g_z. See papers: http://dx.doi.org/10.1364/OE.19.025881
%     and http://dx.doi.org/10.1364/OE.19.012066.
%   'qpcut' : quasiparticle phase retrieval up to a cutoff frequency in
%     Fourier space, see parameter below
%   'qp2' : Similar to the quasiparticle phase retrieval except that the
%     frequency bands are not set to zero but multiplied by the constant
%     which is given by the values of the CTF at the boundary of frequency
%     band to be filtered.
%
% filter_size : 1x2-vector. Default [1024 1024]. Size of output filter.
%
% energy_distance_pixelsize : 1x3-vector. Default: [20e3 0.945 .75e-6].
%   Energy in eV, Distance in m, Pixelsize in m.
%   if imaginary interpreted as Fresnel number.
%
% regularization_parameter : scalar,  default: 2.5. Phase retrieval is
%   regularized according: 
%        1/func(x) -> 1/(func(x)+10^(-regularization_parameter)),
%      'ict': 
%   for details of the placeholder function 'func' see code below. The
%   regularization parameter is the negative of the decadic logartihm of
%   the constant which is added to the denominator in order to regularize
%   the singularity at zero frequency. Typical values are between 1.5 and
%   3.5 depending on energy, residual absorption, etc.
%   Relation to the refractive index is 10^r = delta / beta
%
% binary_filter_threshold : scalar in [0 1], default: 0.1. Parameter for
%   quasiparticle % phase retrieval defining the width of the rings which
%   are cropped around the zero crossings of the CTF denominator (in
%   Fourier space). Typical values are between 0.01 and 0.1, where large
%   values yields results similiar to 'tie' phase retrieval. Note that a
%   value of e.g. 0.1 filters 10% of all frequencies in Fourier space.
%
% frequency_cutoff : scalar, in rad, default pi, cutoff frequency in radian for
%   'qphalfsine' method 
%
% precision : string. 'single' (default), or 'double'.
% 
% Written by Julian Moosmann.
%
% [fourier_filter, parameter_string] = PhaseFilter(method, filter_size, energy_distance_pixelsize, regularization_parameter, binary_filter_threshold, frequency_cutoff, precision)


%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    method = 'tie';
end
if nargin < 2
    filter_size = [1024 1024];
end
if nargin < 3
    energy_distance_pixelsize = [20e3 0.945 .75e-6];
end
if nargin < 4
    regularization_parameter = 2.5;
end
if nargin < 5
    binary_filter_threshold = 0.1;
end
if nargin < 6 
    frequency_cutoff = 1*pi;
end
if nargin < 7
    precision = 'single';
end

%% TODO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO: frequency cut off for all method with no cut off as default
% TODO: consistence of reg par and duality parameter
% TODO: clean up

%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call 'PhaseFilterDual' if duality is assumed. Then, 'regPar' is
% interpreted as -log10 of duality factor epsilon.
strInd = strfind(lower(method),'dual');
switch ~isempty(strInd)
%% Standard phase retrieval %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 0
%% Parameter
if isreal(energy_distance_pixelsize)
    Energy    = energy_distance_pixelsize(1);
    Distance  = energy_distance_pixelsize(2);
    Pixelsize = energy_distance_pixelsize(3);
    lambda    = 6.62606896e-34*299792458/(Energy/1000*1.60217733e-16);
    ArgPrefac = pi*lambda*Distance/Pixelsize^2;

else
    ArgPrefac = pi / abs(energy_distance_pixelsize);
end

%% Fourier coordinates
% 1D
xi  = FrequencyVector(filter_size(2),precision,1);
eta = FrequencyVector(filter_size(1),precision,1);
% 2D
[xi2, eta2]   = meshgrid(xi,eta);
% Function on 2D
sinArg = ArgPrefac*(xi2.^2 + eta2.^2);

%% Filter 
switch lower(method)
    case 'tie'
        fourier_filter = 1/2*sign(sinArg)./(abs(sinArg)+10^-regularization_parameter);
        if regularization_parameter < 0
            parameter_string = sprintf('tie_regPar_m%05.2f',abs(regularization_parameter));
        else
            parameter_string = sprintf('tie_regPar_p%05.2f',regularization_parameter);
        end
    case 'ict'
        sinxiquad   = cos(sinArg) + 10^regularization_parameter * sin(sinArg);
        fourier_filter = 1/2*sign(sinxiquad)./(abs(sinxiquad)+10^-regularization_parameter);        
        parameter_string = sprintf('ict_regPar%3.2f',regularization_parameter);
    case 'ctf'
        sinxiquad   = sin(sinArg);
        fourier_filter = 1/2*sign(sinxiquad)./(abs(sinxiquad)+10^-regularization_parameter);        
        parameter_string = sprintf('ctf_regPar%3.2f',regularization_parameter);
    case {'ctfhalfsine','ctffirsthalfsine','halfsine','firsthalfsine'}
        sinxiquad   = sin(sinArg);
        fourier_filter = 1/2*sign(sinxiquad)./(abs(sinxiquad)+10^-regularization_parameter);
        fourier_filter( sinArg >= pi ) = 0;
        parameter_string = sprintf('ctfHalfSine_regPar%3.2f',regularization_parameter);
    case {'qp','pctf','quasi','quasiparticle','quasiparticles'}
        sinxiquad   = sin(sinArg);
        fourier_filter = 1/2*sign(sinxiquad)./(abs(sinxiquad)+10^-regularization_parameter);
        fourier_filter( sinArg > pi/2  &  abs(sinxiquad) < binary_filter_threshold) = 0;
        parameter_string = sprintf('qp_regPar%3.2f_binFilt%3.3f',regularization_parameter,binary_filter_threshold);
    case {'qpcut', 'quasicut', 'quasiparticlecut'}
        sinxiquad   = sin(sinArg);
        fourier_filter = 1/2*sign(sinxiquad)./(abs(sinxiquad)+10^-regularization_parameter);
        fourier_filter( sinArg > pi/2  &  abs(sinxiquad) < binary_filter_threshold) = 0;
        fourier_filter( sinArg >= frequency_cutoff ) = 0;
        parameter_string = sprintf('qpcut_regPar%3.2f_binFilt%3.3f_cutoff%3.2fpi',regularization_parameter,binary_filter_threshold, frequency_cutoff/pi);
    case {'qp2','quasi2','pctf2'}
        sinxiquad   = sin(sinArg);
        fourier_filter = 1/2*sign(sinxiquad)./(abs(sinxiquad)+10^-regularization_parameter);
        mask = sinArg > pi/2  &  abs(sinxiquad) < binary_filter_threshold;
        fourier_filter( mask ) = bsxfun(@(a,b) a(b),sign(fourier_filter)/(2*(binary_filter_threshold+10^-regularization_parameter)),mask);
        parameter_string = sprintf('qp2_regPar%3.2f_binFilt%3.3f',regularization_parameter,binary_filter_threshold);
end

% Replace dots by p
parameter_string = regexprep(parameter_string,'\.','p');

%% Dual phase retrieval %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 1
        % remove string 'dual' from Variable 'method' to use it as input
        % for 'PhaseFilterDual'
        method(strInd+(0:3)) = [];
        % Call 'PhaseFilterDual'
        [fourier_filter,parameter_string] = PhaseFilterDual(method,filter_size,energy_distance_pixelsize,10^-regularization_parameter,binary_filter_threshold,precision);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function k = FrequencyVector(len, precision, normalized)
% (Normalized) frequency vector matching MATLAB's discretization of Fourier
% space.
%
% ARGUMENTS:
% len: scalar integer. Length of the frequency vector.
% precsision: single' or 'double'. Default: 'single'
% normalized: boolean, default: 1. Output frequencies are normalized
%   between [-1 1].
%
% OUTPUT:
% k : 1D frequency vector, to be used with meshgrid for example.
%
% Written by Julian Moosmann.
%
% k = FrequencyVector(len, precision, normalized)

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    precision = 'single';
end
if nargin < 3
    normalized = 1;
end

%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% normalized
len = eval( sprintf('%s(%u)', precision, len) );

% Create vector
%k = [0:1:ceil( (len-1)/2 )  -floor( (len-1)/2):1:-1] /(1 -normalized*(1 - 1*len ) );
if mod( len, 2 ) == 0    
    k = [0:len / 2 - 1, -len / 2:-1];
else
    k = [0:floor( len / 2), -floor( len / 2 ):-1];    
end
if normalized
    k = k / floor( len / 2);
end
