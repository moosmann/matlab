function N = FresnelNumber(energy_keV,distance_m,apertureRadius_m,pixelsize_m,printInfo)
% Compute Fresnel number.
%
% energy_keV: energy of wave in keV
% distance_m: in m. sample-detector distande z
% apertureRadius_m: in m. Radius of aperture of a screen obstructing the
% incident wave
% pixelsize_m: in m. Effective detector pixel size
%
% Written by Julian Moosmann, 2013-10-15


if nargin < 4
    pixelsize_m = 1e-6;
end
if nargin < 5
    printInfo = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda = EnergyConverter(energy_keV);

N = apertureRadius_m^2/(lambda*distance_m);
N2 = pixelsize_m^2/(lambda*distance_m);

if printInfo
    fprintf('\n')
    fprintf(' Fresnel number: N = a^2/(lambda*z) = %g\n',N);
    fprintf(' Fresnel number: N = pixelsize_m^2/(lambda*z) = %g\n',N2);
    fprintf(' Fresnel diffraction: N >~ 1\n Fraunhofer diffraction: N << 1\n Geometrical optics: N >> 1\n');
    fprintf('\n')
    fprintf(' Energy: E = %g keV\n',energy_keV)
    fprintf(' Wave length: lambda = %g m = %g Angstrom\n',lambda,lambda*1e10);
    fprintf(' distance_m: z = %g m = %g mm\n',distance_m,1000*distance_m)
    fprintf(' Aperture radius: a = %g m = %g micron\n',apertureRadius_m,1e6*apertureRadius_m)    
    fprintf(' sqrt(lambda*z) = %g m = %g micron\n',sqrt(lambda*distance_m),1e6*sqrt(lambda*distance_m))
    fprintf(' 4*(pixel size)^2/lambda*z) = %g\n',4*pixelsize_m^2/lambda*distance_m)
    fprintf('\n')
end



