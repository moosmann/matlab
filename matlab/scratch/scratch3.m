a = zeros( [1, 1000]);
offset = 500;
x = offset + (300:650);
x = mod( x, numel( a) ) + 1;
a(x) = 1;

af = fft( a );

af = af / af(1);


af_a = abs( af );
af_p = phase( af );
af_ifc = abs( af(50:100) );

subplot( 3,1,1)
plot( a, '.' )

subplot( 3,1,2)
plot( af_a, '.' )

subplot( 3,1,3)
plot( af_p, '.' )


fprintf( '\n A : %g', sum( af_a ) )
fprintf( '\n A IFC : %g', sum( af_ifc ) )

fprintf( ' \n ' )