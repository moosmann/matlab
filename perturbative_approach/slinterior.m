function [phase,data,phi] = slinterior(dist,slice,zoom,alpha,padding,padvalue,renormalize)

if (nargin < 7) || isempty(renormalize),renormalize = 1;end;
if (nargin < 6) || isempty(padvalue),padvalue = 0;end;
if (nargin < 5) || isempty(padding),padding =1;end;
if (nargin < 4) || isempty(alpha),alpha = 10;end; 
if (nargin < 3) || isempty(zoom),zoom = 0;end;
normpow = 1;

% Flags for pop-up figures.
if 1,
    show_raw_data      = 1;
    show_phase_maps    = 1;
    show_error_maps    = 1;
    show_spectral_line = 256;
end;
    
% Read data
cd ~/data/shepp-logan_interior;
file_name_prefix = sprintf('sl_E14keV_z%05ucm_res512_os1_',dist);
phase_name_postfix = '-usample-arg2pi.gpbin';
[phase,phase,phase] = loadgpbin(sprintf('%s%03u%s',file_name_prefix,slice, ...
                                phase_name_postfix));
phase = (2*pi*phase(:,end:-1:1))/(2*pi);
[header,data] = pmedf_read(sprintf('%s%03u.edf',file_name_prefix,slice));
% Clip raw data and phase to zoom in the inerior of the phantom.
if zoom 
    [dimx,dimy] = size(phase);
    x = ceil(1/2+dimx/2-zoom/2):floor(dimx/2+zoom/2);
    y = x;
    phase = phase(x,y);
    data  = data(x,y);
end;
% Phase retrieval.
phi = rec(1*data,padding,alpha,renormalize);
if renormalize
    fprintf(1,['Reconstructed with alpha=%u, padding=%u and ' ...
               'renormalized.\n'],alpha,padding);
else
    fprintf(1,['Reconstructed with alpha=%u, padding=%u.\n'],alpha, padding);
end;
% Show exact phase map and raw data.
if show_raw_data,
figure('Name','Exact Phase Map and Raw Data'), ...
imshow([normat(data),normat(phase)], ...
       [],'InitialMagnification','fit'),colorbar;
cd ~/data/shepp-logan_interior/pics;
saveas(gcf,sprintf('%sExactPhase_Data.eps',file_name_prefix),'psc2');
end;
% Error maps, mean error and maximum error in real space.
error_bro         = (abs(normat(phase).^normpow-normat(phi(:,:,1)).^normpow));
error_brocor      = (abs(normat(phase).^normpow-normat(phi(:,:,1)+phi(:,:,2)+phi(:,:,3)+phi(:,:,4)).^normpow));
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
% Error maps, mean error and maximum error in Fourier space.
phase_ft             = fftshift(fft2(normat(phase)));
phi0_ft              = fftshift(fft2(normat(phi(:,:,1))));
phi012_ft            = fftshift(fft2(normat(phi(:,:,1)+phi(:,:,2)+phi(:,:,3)+phi(:,:,4))));
error_bro_ft         = abs((phase_ft)-(phi0_ft));
error_brocor_ft      = abs((phase_ft)-(phi012_ft));
max_error_bro_ft     = max(error_bro_ft(:));
max_error_brocor_ft  = max(error_brocor_ft(:));
mean_error_bro_ft    = mean(error_bro_ft(:));
mean_error_brocor_ft = mean(error_brocor_ft(:));
% Print maximum and mean values of error maps.
fprintf(1,['Maximum Fourier space error: BRO=%.8g, BROCOR=%g (%2.3f%%)\n' ...
           'Mean Fourier space error: BRO=%.8g, BROCOR=%g (%2.3f%%)\n'], ...
        max_error_bro_ft,max_error_brocor_ft,100*max_error_brocor_ft/ ...
        max_error_bro_ft, mean_error_bro_ft, mean_error_brocor_ft,100* ...
        mean_error_brocor_ft/ mean_error_bro_ft);
% Show Images:
% Show phase maps.
if show_phase_maps
    figure('Name','Phase maps: top row: EXACT, BROCOR; bottom row: BRO, COR'), ...
    imshow([normat(phase), ...
            normat(phi(:,:,1)+phi(:,:,2)+phi(:,:,3)+phi(:,:,4)); ...
            normat(phi(:,:,1)), ...
            normat(phi(:,:,2)+phi(:,:,3)+phi(:,:,4))], ...
           [],'InitialMagnification','fit'),colorbar;
    saveas(gcf,sprintf('%sExactPhase_BroCor__Bro_Cor.eps',file_name_prefix),'psc2');
end;
% Show real space and Fourier space error maps.
if show_error_maps
    ErrorMaps(error_bro,error_brocor,file_name_prefix);
    %FourierErrorMaps(error_bro_ft,error_brocor_ft);
end;
% Show spectral line cuts.
if show_spectral_line 
    SpectralPlots(phase,phi);
end;