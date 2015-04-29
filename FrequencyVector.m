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
