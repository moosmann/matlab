function [int, intPadded] = Propagation2(PhaseShift,Absorption,EnergyDistancePixelsize,Padding,printInfo) 
% Intensity according Fresnel propagation for a monochromatic plane wave of
% energy 'EnergyDistancePixelsize'(1), a sample to detector distance
% 'EnergyDistancePixelsize'(2), and detector pixelsize of
% 'EnergyDistancePixelsize'(3) using the 2D maps 'PhaseShift'
% and 'Absorption' as input. Here 'Absorption' is two times the linear
% attenuation, such that the input wave I(z=0) at zero distance is given by
% I(z=0) = exp( i*'PhaseShift' - 'Absorption' ).
%
%
% Written by Julian Moosmann, last version 2013-11-14
%
% [int, intPadded] = Propagation2(PhaseShift,Absorption,EnergyDistancePixelsize,Padding,printInfo) 

%% Defaults %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 3
    EnergyDistancePixelsize = [20 0.945 .75e-6];% [keV m m]
end;
if nargin < 4
    Padding = '';'symmetric';
end;
if nargin < 5
    printInfo = 0;
end;
%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
precision = class(PhaseShift);
normalize = 1;
Energy    = EnergyDistancePixelsize(1);
Distance  = EnergyDistancePixelsize(2);
Pixelsize = EnergyDistancePixelsize(3);
lambda    = 6.62606896e-34*299792458/(Energy*1.60217733e-16);
% Prefactor needed for TIE and CTF retrieval.
ArgPrefac = 2*pi*lambda*Distance/Pixelsize^2;
[dimx, dimy] = size(PhaseShift);

% Pad object
if Padding
    padX = ceil(dimx/2);
    padY = ceil(dimy/2);    
    if ~isscalar(PhaseShift)
        PhaseShift  = padarray(PhaseShift,[padX padY],Padding,'both');
    end
    if ~isscalar(Absorption)
        Absorption  = padarray(Absorption,[padX padY],Padding,'both');
    end
    [dimx, dimy] = size(PhaseShift);
end

% Fourier cooridnates
[xi,eta] = meshgrid(FrequencyVector(dimy,precision,normalize),FrequencyVector(dimx,precision,normalize));

% Intensity according Fresnel propagation
if isscalar(Absorption) && Absorption > 0
    int = abs( ifft2( exp( -1i* ArgPrefac * (xi.^2 + eta.^2)/2 ) .* fft2(exp( (-Absorption + 1i) * PhaseShift)) ) ).^2;    
else
    int = abs( ifft2( exp( -1i* ArgPrefac * (xi.^2 + eta.^2)/2 ) .* fft2(exp( -Absorption + 1i*PhaseShift)) ) ).^2;
end

% Optional output of padded and not yet cropped intensity
if nargout > 1
    intPadded = int;
end

% Crop to orginial size
if Padding
    int = int(padX+1:end-padX,padY+1:end-padY);
end

% Print info
if printInfo,
    fprintf(' \n FRESNEL PROPAGATION: \n')
    fprintf(' E = %g keV, z = %g m, dx = %g micron, size = %u x %u, 2*pi*lambda*z/dx^2 = %g\n', ...
        Energy,Distance,Pixelsize,dimx,dimy,ArgPrefac);
    if Padding
        fprintf(' Padding: %s, Padded size: %u x %u \n',Padding,dimx,dimy)
    end
    domain(PhaseShift)
    if ~isscalar(Absorption)
        domain(Absorption)
    end
    domain(int,1,'Intensity')
    fprintf(' \n')
end;
