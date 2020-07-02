function im = FilterHisto( im, number_of_stds, roi)
% Filter histogram of image such that all values lie within +/-
% 'number_of_stds' standard deviations around the image mean.
% 
% Arguments
% im : 2D matrix
% number_of_std : scalar, default: 3
% roi : scalar in [0 1]. defines region of the image which is used for the
% histogram.
%
% Wirtten by Julian Moosmann. First version: 2017-02-08, last: 2017-02-08

%% Default arguments
if nargin < 2
    number_of_stds = 3;
end
if nargin < 3
    roi = 0;
end

%% Main 
if number_of_stds == 0
    return
end

[d1, d2] = size( im );
if roi > 0
    x = IndexParameterToRange( [roi, 1-roi], d1 );
    y = IndexParameterToRange( [roi, 1-roi], d2 );
else
    x = 1:d1;
    y = 1:d2;
end

im_mean = mean2( im(x,y) );
im_std = std2( im(x,y) );
im_min = min2( im(x,y) );
im_max = max2( im(x,y) );

hist_min = max( im_mean - number_of_stds * im_std, im_min);
hist_max = min( im_mean + number_of_stds * im_std, im_max);

im( im < hist_min ) = hist_min;
im( im > hist_max ) = hist_max;
