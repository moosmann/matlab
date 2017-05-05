function im = MaskingDisc( im, radial_fraction, value)
% Keep central ellipse within the 'radial_fraction' of the largest ellipse
% fitting within the rectangular image and set the region outside of this
% ellipse to the mean within the disc or to 'value'.
%
% im: 2D-matrix
% radial_fraction: scalar in [0 1]. default: 0.95. defines the radius of the
% disc to be keep as a fraction of the largest ellipsis that fits within the
% rectangular image
% value : scalar or []. default: []. if [] uses the mean of disc
%
% Written by Julian Moosmann, last version:2017-04-05
%
% im = MaskingDisc( im, radial_fraction, value) 

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    im = rand( 10, 14 );
end
if nargin < 2
    radial_fraction = 0.95;
end
if nargin < 3
    value = [];
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x, y] = meshgrid(-size(im,2)/2:size(im,2)/2-1,-size(im,1)/2:size(im,1)/2-1);

x = x - (x(1,1)+x(1,end))/2;
y = y - (y(1,1)+y(end,1))/2;

m = sqrt((x/x(1,end)).^2 + (y/y(end,1)).^2);

m = m < radial_fraction;

if isempty(value)
    % set value to inner mean
    value = sum(im(:).*m(:))/sum(m(:));
end

im = m.*im + (1-m)*value;