function test_sliced_variable()

phase_retrieval.method = 'tie';
phase_retrieval.use_parpool = 0;

fprintf( '\n Create array' )
proj = ones( [2001,2000, 2000], 'single');
fprintf( '\n size: %u, %u, %u, %u, %u', size(proj));
proj(1,1) = 10;
proj(end,end) = 0;
fprintf( '\n Call myfunc' )
proj = test_sliced_variable_myfunc(proj, phase_retrieval);
fprintf( '\n size: %u, %u, %u, %u, %u', size(proj));

% fprintf( '\n Create array' )
% a = ones( [4001,4000, 2000], 'single');
% a(1,1) = 10;
% a(end,end) = 0;
% fprintf( '\n Call myfunc' )
% [a, s1, s2] = test_sliced_variable_myfunc(a, 1, 2);
% fprintf( '\n size: %u, %u, %u, %u, %u', size(a), s1, s2 );

end




