function [phase,data,phi0,phi1,phi2] = sl(dist,slice,zoom,alpha,padding,padvalue,renormalize)

if (nargin < 7) || isempty(renormalize), renormalize = 1;
end;
if (nargin < 6) || isempty(padvalue),  padvalue = 0;
end;
if (nargin < 5) || isempty(padding),  padding =1;
end;
if (nargin < 4) || isempty(alpha),  alpha = 10;
end; 
if (nargin < 3) || isempty(zoom),  zoom = 0;
end;
% Flags for pop-up figures.
if 1
    show_raw_data      = 1;
    show_phase_maps    = 1;
    show_error_maps    = 1;
    show_spectral_line = 256;
else
    show_raw_data      = 0;
    show_phase_maps    = 1;
    show_error_maps    = 0;
    show_spectral_line = 0;
end;
% Read data
cd ~/data/sl7/feb03_2330_phaseshiftlesspi;
file_name_prefix = sprintf('sl_E14keV_z%05ucm_res512_os1_',dist);
phase_name_postfix = '-usample-arg2pi.gpbin';
[x,y,phase] = loadgpbin(sprintf('%s%03u%s',file_name_prefix,slice, ...
                                phase_name_postfix));
phase = phase(:,end:-1:1);
data = pmedfread(sprintf('%s%03u.edf',file_name_prefix,slice));
% Show exact phase map and raw data.
if show_raw_data
figure('Name','Exact Phase Map and Raw Data'), ...
imshow([normat(data),normat(phase)], ...
       [],'InitialMagnification','fit'),colorbar;
end;
% Clip raw data and phase to zoom in the inerior of the phantom.
if zoom 
    [dimx,dimy] = size(phase);
    x = ceil(1/2+dimx/2-zoom/2):floor(dimx/2+zoom/2);
    y = x;
    phase = phase(x,y);
    data  = data(x,y);
end;
% Phase retrieval.
[phi0,phi1,phi2] = rec(1*data,padding,alpha,renormalize);
if renormalize
    fprintf(1,['Reconstructed with alpha=%u, padding=%u and ' ...
               'renormalized.\n'],alpha,padding);
else
    fprintf(1,['Reconstructed with alpha=%u, padding=%u.\n'],alpha, padding);
end;
% Show Images:
% Show phase maps.
if show_phase_maps
    figure('Name','Phase maps: phase_exact, phase_rec_brocor; phase_rec_bro, phase_rec_cor'), ...
imshow([normat(phase),normat(phi0+phi1+phi2);normat(phi0),normat(phi1+phi2)], ...
       [],'InitialMagnification','fit'),colorbar;
end;
% Show error maps.
if show_error_maps
    ErrorMaps(phase,phi0,phi1,phi2);
end;
% Show spectral line cuts.
if show_spectral_line 
    SpectralPlots(phase,phi0,phi1,phi2);
end;