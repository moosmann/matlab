function im = FilterHisto( im, number_of_stds)
% Filter histogram of image such that all values lie within +/-
% 'number_of_stds' standard deviations around the image mean.
% 
% Arguments
% im : 2D matrix
% number_of_std : scalar, default: 3
%
% Wirtten by Julian Moosmann. First version: 2017-02-08, last: 2017-02-08

%% Default arguments
if nargin < 2
    number_of_stds = 3;
end

%% Main 
im_mean = mean2( im );
im_std = std2( im );
im_min = min2( im );
im_max = max2( im );

hist_min = max( im_mean - number_of_stds * im_std, im_min);
hist_max = min( im_mean + number_of_stds * im_std, im_max);

im( im < hist_min ) = hist_min;
im( im > hist_max ) = hist_max;