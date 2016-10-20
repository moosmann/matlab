function array = NegLog( array, toggle)
% Take the negative logarithm of the input array.

%% Defaults
if nargin < 2
    toggle = 1;
end

%% Main
if toggle
    array = - reallog( array );
else
    return 
end