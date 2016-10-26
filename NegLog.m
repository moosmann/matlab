function array = NegLog( array, toggle, enforce_positivity)
% Take the negative logarithm of the input array.

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