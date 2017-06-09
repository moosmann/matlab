% Tomographic reconstruction of Shepp-Logan phantom with filtered
% backprojection using MATLAB's iradon and the ASTRA toolbox

pixelsize = 1e-6; %m

%% Phantom
N = 128;
num_proj = 2*N;
theta = 0:180/num_proj:180-1/num_proj;
P = 255e6*phantom('Modified Shepp-Logan',N) * pixelsize;

%% Projections
[R, Xp] = radon(P,theta);
RotAxis = ceil(size(R,1)/2);

%% Inverse radon transformation
freqScal = 1;
OutputSize = size(R,1);
[I,h] = iradon(R,theta,'linear','Ram-Lak',freqScal,OutputSize);


%% ASTRA

% Ramp filter
filt = iradonDesignFilter('Ram-Lak', size(R,1), 1);

% Butterworth filter
[b, a] = butter(1, 0.5);
bw = freqz(b, a, numel(filt) );
%filt = filt .* bw;

% Apply filters
Rf = real( ifft( bsxfun(@times, fft( R, [], 1), filt), [], 1, 'symmetric') );
A = astra_parallel3D( Rf, pi/2 + theta * pi/180, 0, [OutputSize OutputSize 1]);


%% Print %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( '\n shape projections : %u %u', size( R ) )
fprintf( '\n shape phantom : %u %u', size( P ) )
fprintf( '\n shape reco iradon: %u %u', size( I ) )
fprintf( '\n shape reco astra : %u %u', size( A ) )

fprintf( '\n\n DATA RANGE\n' )
domain( R , 1, 'projections')
domain( P , 1, 'phantom    ')
domain( I , 1, 'reco iradon')
domain( A , 1, 'reco astra ')

% 
offset = round( (size(R,1) - N) / 2 );
xx = offset + (1:N);
yy = xx;

%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y1=R(xx,1);
y1=y1(:); 
y2 = sum(P,1);
y2=y2(:);
Y = [y1, y2];

m = 2;n = 3;
nn = 1;
%subplot(m,n,1), imsc(P), title('phantom'), axis image

subplot(m,n,nn), imsc(R), title('projections'), axis image, nn=nn+1;
subplot(m,n,nn), imsc(I(xx,yy)), title('reco iradon'), axis image, nn=nn+1;
subplot(m,n,nn), imsc(I(xx,yy) - P), title('reco iradon - phantom'), axis image, nn=nn+1;colorbar


subplot(m,n,nn), plot( Y ), title('iradon and sum'), legend('radon', 'sum'), nn=nn+1;
subplot(m,n,nn), imsc(A(xx,yy)), title('reco astra'), axis image, nn=nn+1;
subplot(m,n,nn), imsc(A(xx,yy) - P), title('reco astra - phantom'), axis image, nn=nn+1;colorbar

%% adjoint %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf( '\n ADJOINT radon/iradon' )
x = ones( N, N);
Ax = radon( x, theta );
y = ones( size(Ax) );
Ady  = iradon(y,theta,'linear','none',1,N);
Ady = Ady * num_proj/pi * 2;
n1 = sum(Ax(:) .* y(:));
n2 = sum(Ady(:) .* x(:));
fprintf('\n <A x,y>  = %g',n1)
fprintf('\n <x,Ad y> = %g',n2)
fprintf('\n <A x,y> / <x,Ad y> = %g',n1/n2)
fprintf('\n <A x,y> / <x,Ad y> - 1= %g',n1/n2 - 1)

fprintf( '\n\n ADJOINT radon/astra' )
x = ones( N, N);
Ax = radon( x, theta ) ;
y = ones( size(Ax) );
M = 2 * N;
x = ones( M, M);
Ady = astra_parallel3D( y, pi/2 + theta * pi/180, 0, [M M 1], [-N/2 N/2 -N/2 N/2 -0.5 0.5]) ;
Ady = Ady * num_proj/pi * 2;
n1 = sum(Ax(:) .* y(:));
n2 = sum(Ady(:) .* x(:));
fprintf('\n <A x,y>  = %g',n1)
fprintf('\n <x,Ad y> = %g',n2)
fprintf('\n <A x,y> / <x,Ad y> = %g',n1/n2)
fprintf('\n <A x,y> / <x,Ad y> - 1= %g',n1/n2 - 1)

fprintf( '\n\n ADJOINT phantom' )
x = P;
Ax = R;
y = R;
Ady  = I(xx,yy);
Ady = Ady * num_proj/pi * 2;
n1 = sum(Ax(:) .* y(:));
n2 = sum(Ady(:) .* x(:));
fprintf('\n <A x,y>  = %g',n1)
fprintf('\n <x,Ad y> = %g',n2)
fprintf('\n <A x,y> / <x,Ad y> = %g',n1/n2)
fprintf('\n <A x,y> / <x,Ad y> - 1= %g',n1/n2 - 1)



fprintf( '\n' )