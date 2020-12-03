function test_inplace_operation( vol_mem_GB)
%dbstop if error
s = round( (vol_mem_GB*(1024)^3 / 8 )^(1/3) );

fprintf( '\nCreate volume\n ' )
tic
a = 0.1 * ones( [s s s] );
toc
fprintf( ' volume created: %g GB, size : %u %u %u', GB(a), size( a) )


%%
fprintf( '\n\nCalculate mean:\n ' )
tic;
amean = mean( a(:));
toc
fprintf( ' mean : %g', amean )

%%
fprintf( '\n\nStart of: in place multiplication with scalar\n ' )
tic;
a = 3 * a;
toc;
fprintf( ' end of:   in place multiplication with scalar ' )

%%
fprintf( '\n\nStart of: in place multiplication with vector\n ' )
b = 2 * ones( [1 s 1] );
tic;
a = b .* a;
toc;
fprintf( ' end of:   in place multiplication with vector ' )

%
%%
fprintf( '\n\nStart of: in place multiplication with mtimes\n ' )

b = 2 * ones( [1 s 1] );
tic;
a = fun_times( a, b );
toc;
fprintf( ' end of:   in place multiplication with mtimes ' )



%%
fprintf( '\n\nStart of: in place multiplication with parloop\n ' )
b = 2 * ones( [1 s] );
tic;
mem_free = free_memory;
mem_vol = Bytes( a );
num_worker = floor( mem_free / mem_vol );
fprintf( ' num worker max : %u\n', num_worker )
parfor (n = 1:size( a, 3), num_worker )
    im = a(:,:,n);
    im = b .* im + 1;
    a(:,:,n) = im;
end
toc;
fprintf( ' end of:   in place multiplication with parloop ' )

%%
fprintf( '\n\nCalculate mean\n ' )
tic
amean = mean( a(:));
toc
fprintf( ' mean : %g', amean )

fprintf( '\n FIN \n\n' )
end