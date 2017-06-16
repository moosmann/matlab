function mat = normat(mat, MinMax)
% Normalize matrix "mat" to the dynamic range between MinMax(1) and
% MinMax(2), else between [0,1].
%
% Writtten by Julian Moosmann, last version: 2014-06-25

%% Define range
if nargin == 2
    minVal = MinMax(1);
    maxVal = MinMax(2);
else
    minVal = min(mat(:));
    maxVal = max(mat(:));
end

%% Normalize
mat = (mat - minVal)/(maxVal - minVal);
