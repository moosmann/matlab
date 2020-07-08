function A = test_inplace( A, B )

fprintf( '\nbsxfun' )
tic
A = bsxfun( @times, A, B );
toc

fprintf( '\ntimes' )
tic;
A = times( A, B);
toc

fprintf( '\nA.*B' )
tic
A = A .* B;
toc

fprintf( '\nA.*B' )
tic;
A = B .* A;
toc

end