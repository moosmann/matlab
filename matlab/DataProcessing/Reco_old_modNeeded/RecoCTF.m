function [ctf,g_ft,CroppedAreaXY,PaddedDimsXY]=RecoCTF(int,alpha,EnergyDistancePixelsize,BinaryFilter,Padding_FactorAndValue);

% [ctf,g_ft,CroppedAreaXY,PaddedDimsXY]=RecoCTF(int,alpha,EnergyDistancePixelsize,BinaryFilter,Padding_FactorAndValue);
                                                 
if nargin<2,alpha=12;end;
if nargin<3,EnergyDistancePixelsize=[25 1 0.36e-6];end;
if nargin<4,BinaryFilter=0;end;
if nargin<5,Padding_FactorAndValue=[1 0];end;

% PARAMETERS.
PaddingFactor = Padding_FactorAndValue(1);
PaddingValue  = Padding_FactorAndValue(2);
prefactor     = EnergyDistancePixelsize(3)^2/(2*pi*EnergyConverter(EnergyDistancePixelsize(1))*EnergyDistancePixelsize(2));
% Dimensions of input data array.
[dim1,dim2]   = size(int);
xcut          = 1:dim1;
ycut          = 1:dim2;
% Get new dimension whichs are the next-power-of-2 times PaddingFactor. The
% next-power-or-2 is done not to spoil computational power/precision because
% the MATLAB's fft always pads to the next-power-of-2.
dimx         = PaddingFactor*2^nextpow2(dim1);
dimy         = PaddingFactor*2^nextpow2(dim2);
PaddedDimsXY = [dimx dimy];

% PROGRAMME.
% Filters.
[xi,eta]   = meshgrid(-1/2:1/(dimy):1/2-1/(dimy),-1/2:1/(dimx):1/2-1/(dimx));
% xi = matrix of identic row ranging accroding to the first entry of
% meshgrid, the rows are repeated according to the length of the vector
% of the second entry of meshgrid. eta analague to xi.
xi         = fftshift(xi);
eta        = fftshift(eta);
inv_sine   = 1./(2*sign(sin(1/prefactor*(xi.^2 + eta.^2)/2)).*(abs(sin(1/prefactor*(xi.^2 + eta.^2)/2)))+10^-alpha);
% Define Fourier transform of intensity contrast function g=I(x,y,z)/I(x,y,0)-1 and pad it.
% For pure phase objects, <g> states the conservation of flux <I>=1.
int        = padarray(int-mean(int(:)),[ceil((dimx-dim1)/2),ceil((dimy-dim2)/2)],PaddingValue,'pre');
int        = padarray(int-mean(int(:)),[floor((dimx-dim1)/2),floor((dimy-dim2)/2)],PaddingValue,'post');
g_ft       = fft2(int);

%% CTF PHASE RETRIEVAL
% Divide FT[ctf] by the sine in Fourier space, invert
% the Fourier transform, take real part to be on the safe side, substract
% mean since it cannot be retrieved.
ctf       = real(ifft2(inv_sine.*g_ft));
ctf       = ctf-mean(ctf(:));

% Projected CTF
if BinaryFilter>0,
    % Apply binary filter.
    f         = fftshift(SineFilter([dimx dimy],EnergyDistancePixelsize,BinaryFilter));
    % Compute the projected CTF phase.
    pctf      = real(ifft2(inv_sine.*g_ft.*f));
    pctf      = pctf-mean(pctf(:));
    % Crop if padded.
    if dimx~=dim1 | dimy~=dim2,
        xcut  = 1+ceil((dimx-dim1)/2):floor((dimx+dim1)/2);
        ycut  = 1+ceil((dimy-dim2)/2):floor((dimy+dim2)/2);
        ctf   = ctf(xcut,ycut);
        pctf  = pctf(xcut,ycut);
    end;
        ctf   = cat(3,ctf,pctf);
else,
    if dimx~=dim1 | dimy~=dim2,
        xcut  = 1+ceil((dimx-dim1)/2):floor((dimx+dim1)/2);
        ycut  = 1+ceil((dimy-dim2)/2):floor((dimy+dim2)/2);
        ctf   = ctf(xcut,ycut);
    end;
end;
CroppedAreaXY = {xcut; ycut};
