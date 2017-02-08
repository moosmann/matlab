function [fourier_filter, string_appendix] = PhaseFilter(method, filter_size, energy_distance_pixelsize, regularization_parameter, binary_filter_threshold, precision)
% Fourier space filter for Fourier-transform-absed, simple algebraic,
% single-distance phase retrieval.
%
% Given a normalized intensity map I_z (i.e. flat- and dark-field
% corrected), we define the intensity contrast as g_z = I_z - 1. 
%
% Phase retrieval: 
% phi = real( ifft2( fourier_filter(ARG) .* fft2( g_z ) ) ), 
% No fftshift required. Optional output is a unique string consisting of
% the given parameters. 
%
% ARGUMENTS:
%
% method : string. Default: 'tie'.
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
% filter_size : 1x2-vector. Default [1024 1024]. Size of output filter.
%
% energy_distance_pixelsize : 1x3-vector. Default: [20e3 0.945 .75e-6].
% Energy in eV, Distance in m, Pixelsize in m.
%
% regularization_parameter : scalar. Default: 2.5. Phase retrieval is
% regularised according: 1/func(x) ->
% 1/(func(x)+10^(-regularization_parameter)), for details of the
% placeholder  function func see below. Thus, the regularization parameter
% regularization_parameter is the negative of the decadic logartihm of the
% constant which is added to the denominator in order to regularize the
% singularity at zero frequency. Typical values are between 1.5 and 3.5
% depending on energy, residual absorption, etc.
%
% binary_filter_threshold : scalar. Default: 0.1. Parameter for
% quasiparticle % phase retrieval defining the width of the rings which are
% cropped around the zero crossings of the CTF denominator (in Fourier
% space). Typical values are between 0.01 and 0.1, where large values
% yields results similiar to 'tie' phase retrieval.  
%
% precision : string. 'single' (default), or 'double'.
% 
% Written by Julian Moosmann, last modification: 2017-01-06
%
% PhaseFilter(method, filter_size, energy_distance_pixelsize, regularization_parameter, binary_filter_threshold, precision)


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
    precision = 'single';
end

%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call 'PhaseFilterDual' if duality is assumed. Then, 'regPar' is
% interpreted as -log10 of duality factor epsilon.
strInd = strfind(lower(method),'dual');
switch ~isempty(strInd)
%% Standard phase retrieval %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 0
%% Parameter
Energy    = energy_distance_pixelsize(1);
Distance  = energy_distance_pixelsize(2);
Pixelsize = energy_distance_pixelsize(3);
% wave length
lambda    = 6.62606896e-34*299792458/(Energy/1000*1.60217733e-16);
ArgPrefac = pi*lambda*Distance/Pixelsize^2;

%% Fourier coordinates
% 1D
xi  = FrequencyVector(filter_size(2),precision,1);
eta = FrequencyVector(filter_size(1),precision,1);
% 2D
[sinArg, sinxiquad]   = meshgrid(xi,eta);
% Function on 2D
sinArg = ArgPrefac*(sinArg.^2 + sinxiquad.^2);

%% Filter 
switch lower(method)
    case 'tie'
        fourier_filter = 1/2./(sinArg + 10^-regularization_parameter);
        string_appendix = sprintf('tie_regPar%3.2f',regularization_parameter);
    case 'ctf'
        sinxiquad   = sin(sinArg);
        fourier_filter = 1/2*sign(sinxiquad)./(abs(sinxiquad)+10^-regularization_parameter);        
        string_appendix = sprintf('ctf_regPar%3.2f',regularization_parameter);
    case {'ctfhalfsine','ctffirsthalfsine','halfsine','firsthalfsine'}
        sinxiquad   = sin(sinArg);
        fourier_filter = 1/2*sign(sinxiquad)./(abs(sinxiquad)+10^-regularization_parameter);
        fourier_filter( sinArg >= pi ) = 0;
        string_appendix = sprintf('ctfHalfSine_regPar%3.2f',regularization_parameter);
    case {'qp','pctf','quasi','quasiparticle','quasiparticles'}
        sinxiquad   = sin(sinArg);
        fourier_filter = 1/2*sign(sinxiquad)./(abs(sinxiquad)+10^-regularization_parameter);
        fourier_filter( sinArg > pi/2  &  abs(sinxiquad) < binary_filter_threshold) = 0;
        string_appendix = sprintf('qp_regPar%3.2f_binFilt%3.3f',regularization_parameter,binary_filter_threshold);
    case {'qphalfsine','pctfhalfsine','pctffirsthalfsine','quasihalfsine','quasifirsthalfsine'}
        sinxiquad   = sin(sinArg);
        fourier_filter = 1/2*sign(sinxiquad)./(abs(sinxiquad)+10^-regularization_parameter);
        fourier_filter( sinArg > pi/2  &  abs(sinxiquad) < binary_filter_threshold) = 0;
        fourier_filter( sinArg >= pi ) = 0;
        string_appendix = sprintf('qpHalfSine_regPar%3.2f_binFilt%3.3f',regularization_parameter,binary_filter_threshold);
    case {'qp2','quasi2','pctf2'}
        sinxiquad   = sin(sinArg);
        fourier_filter = 1/2*sign(sinxiquad)./(abs(sinxiquad)+10^-regularization_parameter);
        mask = sinArg > pi/2  &  abs(sinxiquad) < binary_filter_threshold;
        fourier_filter( mask ) = bsxfun(@(a,b) a(b),sign(fourier_filter)/(2*(binary_filter_threshold+10^-regularization_parameter)),mask);
        string_appendix = sprintf('qp2_regPar%3.2f_binFilt%3.3f',regularization_parameter,binary_filter_threshold);
end

% Restore zero frequency component
fourier_filter(1) = 1/2*10^regularization_parameter;

% Replace dots by p
string_appendix = regexprep(string_appendix,'\.','p');

%% Dual phase retrieval %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 1
        % remove string 'dual' from Variable 'method' to use it as input
        % for 'PhaseFilterDual'
        method(strInd+(0:3)) = [];
        % Call 'PhaseFilterDual'
        [fourier_filter,string_appendix] = PhaseFilterDual(method,filter_size,energy_distance_pixelsize,10^-regularization_parameter,binary_filter_threshold,precision);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function k = FrequencyVector(len, precision, normalized)
% Frequency vector of length N corresponding to MATLAB's implementation of
% disrectized Fourier space discretization.
%
% len: scalar integer
% precsision: single' or 'double'. Default: 'single'
% normalized: boolean, default: 1. Output frequencies are normalized between
% [-1 1].
%
% Written by Julian Moosmann, last version 2013-11-12. Update: 2017-10-26
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

% function k = FrequencyVector(SizeOfVector,Precision,Normalize)
% % Frequency vector for MATLAB's implementation of Fourier space
% % discretization (without 'fftshift') ranging from 0 to
% % floor(('SizeOfVector'-1)/2) and -floor('SizeOfVector'/2) to -1.
% %
% % SizeOfVector: integer.
% %
% % Precsision: string, default: 'single'. Available: 'single', or 'double'.
% %
% % Normalize: boolean, default: 1. Output frequencies normalized by size of
% % vector.
% %
% % Written by Julian Moosmann, last version 2013-11-12.
% %
% %k = FrequencyVector(SizeOfVector,Precision,Normalize)
% 
% %% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if nargin < 2
%     Precision = 'single';
% end
% if nargin < 3
%     Normalize = 1;
% end
% 
% %% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Precision
% SizeOfVector = eval(sprintf('%s(%u)',Precision,SizeOfVector));
% 
% % Create vector
% k = [0:1:ceil( (SizeOfVector-1)/2 )  -floor( (SizeOfVector-1)/2):1:-1] /(1 -Normalize*(1 - 1*SizeOfVector ) );