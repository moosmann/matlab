function sinoc = CropShift( sino, shift )
% Crop sinograms with lateral shifts such that the rotation axis is
% centered.
%
% Writen by J. Moomann, 2018-08-04. Last version: 2018-08-10

if (numel( shift ) == 1 && shift == 0) || isempty( shift )
    sinoc = sino;
    return
end

x0 = round( 1 + shift' - min( shift(:) ) );
x1 = size(sino,1) - max( x0 ) + x0;

sinoc = zeros( [size(sino,1) - max(x0) + 1, size(sino,2)] ,'single' );

% parfor nn = 1:size( sino, 2)
%     l = sino(:,nn);
%     sinoc(:,nn) = l(x0(nn):x1(nn));
% end

for nn = 1:size( sino, 2)
    sinoc(:,nn) = sino(x0(nn):x1(nn),nn);
end
