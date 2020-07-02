function FourierErrorMaps(error_bro,error_brocor,zoom)
%Show table of images of the differences of the exact phase and the
%reconstructed phases: error_BRO, error_BROCOR

if (nargin<3) || isempty(zoom)
  zoom = 0;
end;
% Clip raw data and phase to zoom in the inerior of the phantom.
if zoom 
    [dimx,dimy] = size(error_bro);
    x = ceil(1/2+dimx/2-zoom/2):floor(dimx/2+zoom/2);
    y = x;
    error_bro = error_bro(x,y);
    error_brocor = error_brocor(x,y);
end;
% Plot error map.
figure('Name','Log maps of phase error in Fourier space (2D): left=EXACT-BRO, right=EXACT-BROCOR'), ...
imshow([log10(error_bro),log10(error_brocor)],[],'InitialMagnification','fit'),colorbar;




