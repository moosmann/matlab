function prefac(energy,dist,pixelsize)

if nargin < 2 && numel(energy) == 3
    dist = energy(2);
    pixelsize = energy(3);
    energy = energy(1);
end

lambda = EnergyConverter(energy/1000);
pf = 2*pi*lambda*dist/(pixelsize)^2;
fprintf('Energy: E = %g keV\n', energy/1000)
fprintf('Wave length: lambda = %g angstrom\n', lambda * 1e10) 
fprintf('Distance: z = %g mm\n', 1000*dist)
fprintf('Pixel size: dx = %.3g micron\n', 1e6 * pixelsize)
fprintf('lambda*z = %8.3g m^2\n',lambda*dist)
fprintf('Argument of the ctf sine: prefactor*(xi^2 + eta^2)/2\n')
fprintf('Range of coordinates in Fourier space: [xi eta] = [-1/2 1/2]\n')
fprintf('prefactor = 2*pi*lambda*z/pixelsize^2 = %06.3f\n',pf)
fprintf('Maximum argument of the sine at [|xi| |eta|] = [1/2 1/2]: %5.3f\n',pf/4)
fprintf('Zero crossings - 1 (without the central one at xi = eta = 0): %5.3f\n',pf/pi/4);
