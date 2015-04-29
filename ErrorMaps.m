function ErrorMaps(error_bro,error_brocor,file_name_prefix,zoom)
%Show table of images of the differences of the exact phase and the
%reconstructed phases: error_BRO, error_BROCOR

if (nargin<3)||isempty(file_name_prefix),file_name_prefix='Data_';end;
if (nargin<4)||isempty(zoom),zoom = 0;end;

% Clip raw data and phase to zoom in the inerior of the phantom.
if zoom 
    [dimx,dimy] = size(error_bro);
    x = ceil(1/2+dimx/2-zoom/2):floor(dimx/2+zoom/2);
    y = x;
    error_bro = error_bro(x,y);
    error_brocor = error_brocor(x,y);
end;
% Plot error map.
figure('Name','Error maps of phase in real space (2D): left=EXACT-BRO, right=EXACT-BROCOR'), ...
imshow([error_bro,error_brocor],[],'InitialMagnification','fit'), ...
    colorbar;
saveas(gcf,sprintf('%sRealSpaceErropMaps__Bro_BroCor.eps',file_name_prefix),'psc2');



