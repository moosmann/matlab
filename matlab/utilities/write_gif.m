function write_gif( vol, filename )

% Normalize slicewise to [0, 1]
slice_min = min( min( vol ) );
slice_max = max( max( vol ) );
sequ = bsxfun( @rdivide, bsxfun( @minus, vol, slice_min ), slice_max - slice_min );

% 8-bit conversion
sequ = uint8( 2^8 * sequ );

% Insert required singleton dimension
sequ = permute( sequ, [1 2 4 3] );

% Save image sequence as gif
imwrite( sequ, filename )