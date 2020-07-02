function array = NegLog( array, toggle, enforce_positivity)
% Take the negative logarithm of the input array.
%
% ARGUMENTS
% array : N-D array
% toggle : bool. default: 1. take negative logarithm or not
% enforce_positivity : bool. default: 0. if 1: negative values are set to 0
%
% Written by Julian Moosmann, last mod: 2017-06-03
%
% array = NegLog( array, toggle, enforce_positivity)

%% Defaults
if nargin < 2
    toggle = 1;
end
if nargin < 3
    enforce_positivity = 0;
end

%% Main
if toggle
    array = - reallog( array );
    if enforce_positivity
        array( array < 0 ) = 0;
    end
else
    return 
end
