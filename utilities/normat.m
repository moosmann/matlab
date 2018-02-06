function array = normat( array, min_max)
% Normalize matrix "array" to the dynamic range between min_max(1) and
% min_max(2), else between [0,1].
%
% Writtten by Julian Moosmann, last version: 2014-06-25

%% Define range
if nargin == 2
    min_val = min_max(1);
    max_val = min_max(2);
else
    min_val = min(array(:));
    max_val = max(array(:));
end

%% Normalize
array = ( array - min_val )/ ( max_val - min_val );
