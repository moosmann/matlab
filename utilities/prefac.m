function prefac(energy,dist,pixelsize)
    

lambda = EnergyConverter(energy);
pf = 2*pi*lambda*dist/(pixelsize)^2;
fprintf('Energy:      E = %.3u keV\n',energy)
fprintf('Wave length: lambda = %.3g m\n',lambda) 
fprintf('Distance:    z = %4.3f m\n',dist)
fprintf('Pixel size:  dx = %.3g m\n',pixelsize)
fprintf('lambda*z = %8.3g m^2\n',lambda*dist)
fprintf('Argument of the ctf sine: prefactor*(xi^2 + eta^2)/2\n')
fprintf('Range of coordinates in Fourier space: [xi eta] = [-1/2 1/2]\n')
fprintf('prefactor = 2*pi*lambda*z/pixelsize^2 = %06.3f\n',pf)
fprintf('Maximum argument of the sine at [|xi| |eta|] = [1/2 1/2]: %5.3f\n',pf/4)
fprintf('Zero crossings - 1 (without the central one at xi = eta = 0): %5.3f\n',pf/pi/4);
