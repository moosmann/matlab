function zercro = prefac(energy_eV, dist_m, pixelsize_m, verbose)

%% Default arguments
if nargin < 2 && numel(energy_eV) == 3
    dist_m = energy_eV(2);
    pixelsize_m = energy_eV(3);
    energy_eV = energy_eV(1);
end
if nargin < 4
    verbose = 1;
end

%% Main
lambda = EnergyConverter(energy_eV/1000);
pf = 2*pi*lambda*dist_m/(pixelsize_m)^2;
zercro = pf/pi/4;

% Print
if verbose
    fprintf('energy_eV: E = %g keV\n', energy_eV/1000)
    fprintf('Wave length: lambda = %g angstrom\n', lambda * 1e10) 
    fprintf('Distance: z = %g mm\n', 1000*dist_m)
    fprintf('Pixel size: dx = %.3g micron\n', 1e6 * pixelsize_m)
    fprintf('lambda*z = %8.3g m^2\n',lambda*dist_m)
    fprintf('Argument of the ctf sine: prefactor*(xi^2 + eta^2)/2\n')
    fprintf('Range of coordinates in Fourier space: [xi eta] = [-1/2 1/2]\n')
    fprintf('prefactor = 2*pi*lambda*z/pixelsize_m^2 = %06.3f\n',pf)
    fprintf('Maximum argument of the sine at [|xi| |eta|] = [1/2 1/2]: %5.3f\n',pf/4)
    fprintf('Zero crossings - 1 (without the central one at xi = eta = 0): %5.3f\n',zercro);
end
