% Plot phase and modulus of a theta function

for offset = [0, 200, 500, 750]
    for width = 250
        
        a = zeros( [1, 1000]);
        x = offset + (1:width);
        x = mod( x, numel( a) ) + 1;
        a(x) = 1;
        
        af = fft( a );
        
        af = af / af(1);
        
        
        af_a = abs( af );
        af_p = phase( af );
        af_ifc = abs( af(50:100) );
        
        figure( 'Name', sprintf( 'offset: %u, width: %u', offset, width ))
        subplot( 3,1,1)
        plot( a, '.' )
        title( 'signal' )
        
        
        subplot( 3,1,2)
        plot( af_a, '.' )
        title( 'abs( fft( signal ) )')
        
        subplot( 3,1,3)
        plot( af_p, '.' )
        title( 'phase( fft( signal ) )')
        
        fprintf( '\n A : %g', sum( af_a ) )
        fprintf( '\n A IFC : %g', sum( af_ifc ) )
        
        fprintf( ' \n ' )
    end
end
