function ErrorPlots2DZoom(phase_exact,phase_rec_0,phase_rec_1,phase_rec_2,show_unzoomed)
%Show table of images of the differences of the exact phase and the
%reconstructed phases (phase_rec_0, phase_rec_1, phase_rec_2)
    
if (nargin < 5) || isempty(show_unzoomed)
  show_unzoomed=0;
end;
%zoom range
[dimx,dimy] = size(phase_exact);
a = dimx/5;
zoom_range_x = ceil(dimx/2 - a):floor(dimx/2 + a); 
zoom_range_y = ceil(dimy/2 - a):floor(dimy/2 + a);
%error maps, unzoomed and zoomed
error_0   = abs(normat(phase_exact)-normat(phase_rec_0));
error_012 = abs(normat(phase_exact)-normat(phase_rec_0+phase_rec_1+ ...
                                           phase_rec_2));
error_zoomed_0   = error_0(zoom_range_x,zoom_range_y);
error_zoomed_012 = error_012(zoom_range_x,zoom_range_y);

%error plots unzoomed
if show_unzoomed
figure('Name','ErrorPlots2D: phase_exact-phase_rec_0,phase_exact-phase_rec_(0+1+2)'), ...
imshow([error_0,error_012],[],'InitialMagnification','fit'),colorbar;
end;
%error plots zoomed
figure('Name','ErrorPlots2D_zoomed: phase_exact-phase_rec_0,phase_exact-phase_rec_(0+1+2)'), ...
imshow([error_zoomed_0,error_zoomed_012], ...
       [],'InitialMagnification','fit'),colorbar;
%print maximum values of error maps
fprintf(1,['maximum error_0: %g and maximum error_012: %g\n' ...
        'maximum error_zoomed_0: %g and maximum error_zoomed_012: %g\n'], ... 
        max(error_0(:)),max(error_012(:)),max(error_zoomed_0(:)),max(error_zoomed_012(:)));



