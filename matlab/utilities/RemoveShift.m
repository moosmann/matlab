function Vq = RemoveShift( sino, offset_shift )
% Remove lateral offset shift by interpolation. First dimension: pixel,
% second dimension: angle.
%
% Written by J. Moosmann.

[num_pix, num_proj] = size( sino );

% Query grid
x_left = ceil( max( offset_shift) ) + 1 ;
x_right =  floor( num_pix + min( offset_shift ) ) - 0 ;
x = x_left:x_right;
y = 1:num_proj;
[Yq, Xq] = meshgrid( y, x );
Xq = Xq + offset_shift';
Xq = Xq - min( Xq(:) ) + 1;


% Interpolation
V = sino;
method = 'linear';
Vq = interp2( V, Yq, Xq, method );

domain( Xq(:) )
imsc( Vq )
m = isnan( Vq(:));
if sum( m ) > 0
    warning( 'NaN after shift removal' )
end

