ca;
clear;

bin = 4;
N0 = 3056;
energy = 12.25e3; 10e3; % eV
z = 1.25;
pixelsize0 =  1 / bin * 5.5e-9; 0.3e-6; %m

N = round( N0 /bin );

pixelsize = bin * pixelsize0;

proj = 0.5 * normat( phantom( 'shepp-logan', N) );


lambda = E_to_lambda(energy);
k = 2 * pi / lambda;
b = N * pixelsize;
NF = b^2 / lambda / z;
fprintf( '\n number pixels : %u', N );
fprintf( '\n number pixels : %g', pixelsize );
fprintf( '\n characteristic length (detector width) b : %f mm', b * 1000 )
fprintf( '\n energy : %g keV', energy / 1000 );
fprintf( '\n lambda : %f pm', lambda*1e12 )
fprintf( '\n k : %g nm^-1', k*1e-9 )
fprintf( '\n propagation distance : %f m', z );
fprintf( '\n effective pixel size binned : %g micron', pixelsize * 1e6);
fprintf( '\n Fresnel number = b^2 / lambda / z: %g', NF )
fprintf( '\n')

%z0 = 0.1;
%z1 = 10000;
numz = 20;
s(numz).z = 0;
for nn = 1:numz
    %znn = z0 + (z1 - z0) / (numz - 1) * (nn - 1);
    znn = 2^nn * 1e-2;
    s(nn).nn = nn;
    s(nn).z = znn;
    s(nn).NF = b^2 / lambda / znn;
    xmax = znn / ( 2 * pixelsize * k );
    s(nn).xmax = xmax;
    fovscal = xmax / b;
    s(nn).fovscal = fovscal;
    s(nn).int = Propagation( proj, [energy znn pixelsize], 1, 'symmetric', 0);
end

fprintf( '\n')
fprintf( ' %8u', [s.nn] )
fprintf( '\n')
fprintf( ' %8.1f', [s.z] )
fprintf( '\n')
fprintf( ' %8.1f', [s.NF] )
fprintf( '\n')
fprintf( ' %8f', [s.xmax] / b )
fprintf( '\n')
v = cat( 3, s.int );
nimplay( v )