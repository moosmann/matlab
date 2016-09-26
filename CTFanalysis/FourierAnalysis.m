function FourierAnalysis(dist,slice,kzero_rings)
% Compute 2D FT of data. Check for rings of diminuised intensity in
% Fourier space at distance k=sqrt(k_x^2+k_y^2) determined by
% lambda*distance*k=n, where n=1,2,3... .
if (nargin < 3) || isempty(kzero_rings)
    kzero_rings = 1;
end;
    
% Read data
cd ~/data/sl7/feb03_2330_phaseshiftlesspi;
file_name_prefix = sprintf('sl_E14keV_z%05ucm_res512_os1_',dist);
file_name = sprintf('%s%03u.edf',file_name_prefix,slice);
%cd ~/data/sl7/feb03_1415;
%file_name = 'sl_E30keV_z900cm_res512_os1_003.edf'
[~,data] = pmedf_read(file_name);clear header;
% Parameters:
% Resolution
[dimx,dimy] = size(data);
% Detector size
dx = 0.001024;
% Incident wavelength in nm corresponding to 14 keV
lambda = 0.08856030*10^-9;
% Minimal squared radius in Fourier space where intensity diminished ring should be
kzero = (dx/dimx)^2/(dist*10^-2*lambda);
% Fourier trafo of data                                   
fdata = fftshift(fft2(data));
% Create matrix of k_x^2+k_y^2
[k_x,k_y]=meshgrid(-dimx/2:dimx/2-1,-dimy/2:dimy/2-1);
r2 = 1/dimx/dimy*(k_x.^2+k_y.^2);
% Compute rings of diminished intensity for different n
ring = zeros(dimx,dimy);
if kzero_rings
delta_kzero = 0.05*kzero;
fprintf(1,'n*kzero: ');
for n=1:20
    ring(find(n*kzero-delta_kzero<r2 & r2<n*kzero+delta_kzero))=0.5;
    fprintf(1,'%g  ',n*kzero);
end;
fprintf(1,'\n');
end;
% Prepare show image
delta_ft = 5;
absring = ring;
absring(find(abs(fdata)<delta_ft))=1;
figure('Name',file_name),imshow([absring],[], ...
              'XData',[-dimx/2,dimx/2], ... 
              'YData',[-dimy/2,dimy/2]),colorbar;

if 0
argring = ring;
argring(find(mod(angle(fdata),pi)<delta_ft))=1;
realring = ring;
ranges(realring);
realring(find(real(fdata).^2+imag(fdata).^2<delta_ft))=1;
imagring = ring;
imagring(find(imag(fdata)<delta_ft))=1;
figure,imshow([absring,argring;realring,imagring],[]),colorbar;
end;





