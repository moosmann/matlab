%% Parameters
deltaz_m = 6e-6; % micron
energy_keV = 200; % ~ 10^2 kV per electron
lambda_m = EnergyConverter(energy_keV);
pixelsize_m = 0.24e-9; % nm
cSpherAbb_m = 2e-3; % mm
N = 300;
imSize = [N N];

%lambda    = 6.62606896e-34*299792458/(Energy*1.60217733e-16);
prefacSta = pi * lambda_m * deltaz_m / pixelsize_m^2;
prefacSpheAbb = pi * cSpherAbb_m * lambda_m^3 / 2 / pixelsize_m^4;

%% Fourier coordinates
% 1D
outputPrecision = 'double';
xi  = FrequencyVector(imSize(2),outputPrecision,1);
eta = FrequencyVector(imSize(1),outputPrecision,1);
% 2D
[xi, eta]   = meshgrid(xi,eta);
% Function on 2D
xi = (xi.^2 + eta.^2);
arg = prefacSta * xi - prefacSpheAbb * xi.^2;

fprintf('\nPrefactor 1: %g',prefacSta)
fprintf('\nPrefactor 1: %g\n',prefacSpheAbb)
fprintf('\nMax argument: %g\n', max( arg(:) ) )
fprintf('\nMax argumen: %g\n', max( arg(:) ) )

imagesc(fftshift(sin(arg)))
colormap(gray)
