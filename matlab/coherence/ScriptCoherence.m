function Coherence()
%% Parameters

edp_2bm = [30 0.620 2.2e-6];
pixelsize = 0.75e-6; %m
energy = 30;% keV
lambda = EnergyConverter(energy);
z = 0.6;% m
numPix = 1000;
aperture = numPix * pixelsize;
cohLen = 30e-6;
%%  Functions
p = @(x) 4*x^2/lambda;
n = @(a) a^2/z/lambda;
%% Output


fprintf('dx / micron:        %4g %4g %4g \n',1e6*pixelsize,1e6*cohLen,1e6*aperture)
fprintf('z = 4 (dx)^2/lamda: %4g %4g %4g \n',p(pixelsize),p(cohLen),p(aperture))
fprintf('N =   : %4g %4g %4g \n',n(pixelsize),n(cohLen),n(aperture))

fprintf('\nFresnel number for different parameters:\n')
PrintN(20,0.945,0.75e-6)
PrintN(20,0.945,100e-6)
PrintN(20,0.945,2000e-6)
PrintN(30,0.62,2.2e-6)
PrintN(30,0.62,10e-6)
PrintN(30,0.62,2000e-6)

fprintf('\nCoherence length for different parameters:\n')
CoherenceLength(30,10.3e-6,145)
CoherenceLength(20,10.3e-6,145)
CoherenceLength(20,57e-6,145)

CoherenceLength(30,30e-6,45)

function PrintN(energy,distance,aperture)
lambda = EnergyConverter(energy);
fprintf('E: %2g keV (%8g nm), ',energy,1e12*lambda)
fprintf('z: %4g mm, ',distance*1000)
fprintf('aperture: %4g micron, ',aperture*1e6)
fprintf('N: %8g\n',aperture^2/distance/lambda)

function CoherenceLength(Energy,SourceSize,SourceSampleDistance)
lambda = EnergyConverter(Energy);
S = SourceSize;
R = SourceSampleDistance;
da = lambda^2 * R^2 / S;
dx = sqrt(da);
dr = sqrt(da/pi);

fprintf('E: %2g keV (%8g nm), ',Energy,1e12*lambda)
fprintf('Source size: %4g micron, ',S*1e6)
fprintf('Source sample : %4g m, ',R)

fprintf('Coherence length (radius): %8g micron (%8g micron)\n',dx*1e6,dr*1e6)
