clear all;
pardir = '/home/jmoosmann/data/tem/TMV';
filename = sprintf('%s/tiltseries_nonoise_def10.mrc', pardir);
[int0, s, mi, ma, av] = ReadMRC(filename);
% itool(int0)

roi = [0 1];
%roi = [0.5 0.85];

%% Parameters
% Acceleration voltage = 200 keV
energy_keV = 200;
lambda_m = EnergyConverter(energy_keV);
%Defocus: 3.0 micrometer under focus 
defocus_m = 10e-6; 
%Focal length: 2.7 mm
focalLength_m = 2.7e-3;
% Magnification: 25000
magnification = 2500;
%Detector pixel size: 16 micrometer
detectPixel_m = 16e-6;
effectPixel_m = detectPixel_m / magnification;
% Spherical aberration: 2.1 mm
cs_m = 2.1e-3;
% Chromatic aberration: 2.1 mm
cc_m = 2.2e-3; 
imSize = [1024 1024];
% Distance object lens
q_m = focalLength_m * ( magnification + 1 );

a = 23;

%lambda    = 6.62606896e-34*299792458/(Energy*1.60217733e-16);
pfDefocus = pi * lambda_m * (a*defocus_m ) / effectPixel_m^2;
pfSpher = 0* pi * cs_m * lambda_m^3 / 2 / effectPixel_m^4;


%% Fourier coordinates
% 1D
outputPrecision = 'double';
xi  = FrequencyVector(imSize(2),outputPrecision,1);
eta = FrequencyVector(imSize(1),outputPrecision,1);
% 2D
[xi, eta]   = meshgrid(xi,eta);
% Function on 2D
xi = (xi.^2 + eta.^2);
xi = fftshift( xi );
arg = pfDefocus * xi - pfSpher * xi.^2;
%figure(1)
%imsc( 1./( 1 + abs( sin( arg ))))
%figure(2)

fprintf('\nDefocus prefactor: %g',pfDefocus)
fprintf('\nSpherical abberation prefactor: %g\n',pfSpher)
fprintf('Max argument: %g\n', max( arg(:) ) )

%% Phase retrieval filter
regPar = 0;
pf = PhaseFilter('ctfdual', imSize, [energy_keV, a*defocus_m, effectPixel_m], regPar);

im1 = log(1+ abs(fftshift(fft2(int0))));
im2 = fftshift(pf);
y = floor(size(int0, 2)/2);
imsc([normat(im1(:,1:y)) normat(im2(:,end-y:end))])