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