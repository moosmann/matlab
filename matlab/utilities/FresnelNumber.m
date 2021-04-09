function N = FresnelNumber(energy,distance,apertureRadius,pixelsize,printInfo)
% Compute Fresnel number.
%
% energy: energy of wave in keV
% distance: in m. sample-detector distande z
% apertureRadius: in m. Radius of aperture of a screen obstructing the
% incident wave
% pixelsize: in m. Effective detector pixel size
%
% Written by Julian Moosmann, 2013-10-15


if nargin < 4
    pixelsize = 1e-6;
end
if nargin < 5
    printInfo = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda = EnergyConverter(energy);

N = apertureRadius^2/(lambda*distance);
N2 = pixelsize^2/(lambda*distance);

if printInfo
    fprintf('\n')
    fprintf(' Fresnel number: N = a^2/(lambda*z) = %g\n',N);
    fprintf(' Fresnel number: N = pixelsize^2/(lambda*z) = %g\n',N2);
    fprintf(' Fresnel diffraction: N >~ 1\n Fraunhofer diffraction: N << 1\n Geometrical optics: N >> 1\n');
    fprintf('\n')
    fprintf(' Energy: E = %g keV\n',energy)
    fprintf(' Wave length: lambda = %g m = %g Angstrom\n',lambda,lambda*1e10);
    fprintf(' Distance: z = %g m = %g mm\n',distance,1000*distance)
    fprintf(' Aperture radius: a = %g m = %g micron\n',apertureRadius,1e6*apertureRadius)    
    fprintf(' sqrt(lambda*z) = %g m = %g micron\n',sqrt(lambda*distance),1e6*sqrt(lambda*distance))
    fprintf(' 4*(pixel size)^2/lambda*z) = %g\n',4*pixelsize^2/lambda*distance)
    fprintf('\n')
end



