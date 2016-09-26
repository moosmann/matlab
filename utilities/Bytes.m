function out = Bytes(array, PowerOf1024)
% Returns the number of bytes allocated by the variable 'array'.
%
% Paramters
% array: MATLAB workspace variable (not as string)
% PowerOf1024: int or string ('b', 'k', 'm', 'g', see code), unit prefix of
% output
%
% Written by J. Moosmann, last version: 2015-08-25


%% Defaults
if nargin < 2
    PowerOf1024 = 0;
end

%% Main
switch lower(PowerOf1024)
    case {'b','bytes'}
        PowerOf1024 = 0;
    case {'k','kb'}
        PowerOf1024 = 1;
    case {'m','mb'}
        PowerOf1024 = 2;
    case {'g','gb'}
        PowerOf1024 = 3;
end

out = whos('array');
out = out.bytes/1024^PowerOf1024;
