function ErrorMaps2(phi,phase,zoom)
%Show table of images of the differences of the exact phase and the
%reconstructed phases: error_BRO, error_BROCOR

if (nargin<3)||isempty(zoom),zoom = 0;end;
normpow = 1;

% Clip raw data and phase to zoom in the inerior of the phantom.
if zoom 
    [dimx,dimy] = size(error_bro);
    x = ceil(1/2+dimx/2-zoom/2):floor(dimx/2+zoom/2);
    y = x;
    error_bro = error_bro(x,y);
    error_brocor = error_brocor(x,y);
end;
% Error maps, mean error and maximum error in real space.
error_bro         = (abs(normat(phase).^normpow-normat(phi(:,:,1)).^normpow));
error_brocor      = (abs(normat(phase).^normpow-normat(phi(:,:,1)+phi(:,:,2)).^normpow));
max_error_bro     = max(error_bro(:));
max_error_brocor  = max(error_brocor(:));
mean_error_bro    = mean(error_bro(:));
mean_error_brocor = mean(error_brocor(:));
% Print maximum and mean values of error maps.
fprintf(1,['Maximum real space error: BRO=%.8g, BROCOR=%g (%2.3f%%)\n' ...
           'Mean real space error: BRO=%.8g, BROCOR=%g (%2.3f%%)\n'], ...
        max_error_bro,max_error_brocor,100*max_error_brocor/max_error_bro, ...
        mean_error_bro, mean_error_brocor,100*mean_error_brocor/ ...
        mean_error_bro);

% Plot error map.
figure('Name','Error maps of phase in real space (2D): left=EXACT-BRO, right=EXACT-BROCOR'), ...
imshow([error_bro,error_brocor],[],'InitialMagnification','fit'), ...
    colorbar;
%saveas(gcf,sprintf('%sRealSpaceErropMaps__Bro_BroCor.eps',file_name_prefix),'psc2');
