% xeno 4-cell stage, 20keV, TIE, slice 597 of 1494, roi

alphaTIE = alphaCTF + log10(2*pi*lambda*Distance/Pixelsize^2);
% Filter for CTF phase retrieval.
prefactor     = Pixelsize^2/(2*pi*lambda*Distance);
SineXiEta   = sin(1/prefactor*(xi.^2 + eta.^2)/2);
InverseSine = 1./(2*sign(SineXiEta).*(abs(SineXiEta))+10^-alphaCTF);
int       = int-mean(int(:));
int       = fft2(int);
out.ctf       = real(ifft2(InverseSine.*int));
out.ctf       = out.ctf-mean(out.ctf(:));
out.tieLO       = 1./(xi.^2 + eta.^2 + 10^-alphaTIE).*int;

a1 = 5;
v1 = -31.97;
min1 = -34.879;
max1 = 4.254;
mean1 = -17.737;
stddev1 = 12.814;

a2 = 5.32;
v2=-52.633;
min2 = -57.637;
max2 = -15.193;
mean2 = -38.503;
stddev2 = 12.814;


v1/v2