N = 4000;
w = 22;
p = 100;
%x = zeros( N, 1);
ind = 0:N-1;
ind = mod( ind, p);
x = ind;
x(x<w) = 1;
x(x>=w) = 0;

fprintf( '\n width: %u', w )
fprintf( '\n period: %u', p )

[a, l] = xcorr( x );
b = a(N:end) ;
%figure('Name', 'Autocorrlation' )
subplot(4,1,1)
plot( x )
title( sprintf( 'signal. period = %u, width = %u', p, w) )
axis tight

subplot(4,1,2)
plot( a )
title( 'autokorrelation full range' )
axis tight

subplot(4,1,3)
plot( b )
title( 'autokorrelation half range' )
axis tight

subplot(4,1,4)
plot( l )
title( 'lag' )
axis tight

% Width
n = 0;
while b(n+1) > 0
    n = n + 1;
end
aw = n;
fprintf( '\n width from AC: %u', aw )

% Period
b2 = b;
b2(1:aw) = 0;
[~,ap] = max(b2);
ap = ap - 1;
fprintf( '\n period from AC: %u', ap )

fprintf( '\n' )

