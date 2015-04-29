function [m,I,phi] = Analysis_Delta_Dot(distance)
% Analysis of delta-dot-like phase map.

cd ~/data/psf
% Parameters.
dot_magnitude = 2*pi*10/17;
dim = 1024;
%distance = 10;
lambda = 0.5166e-10;
padding = 1;
alpha = 15;
erno = 1;
% Phase Map. Intensity pattern. Reconstruction.
m   = Delta_Dot(dim,1,dot_magnitude);
I   = Propagation(m,distance,lambda,padding);
phi = rec(I,padding,alpha,0);
% Analysis.
error_bro         = (abs(normat(m).^erno-normat(phi(:,:,1)).^erno)).^(1/erno);
error_brocor      = (abs(normat(m).^erno-normat(phi(:,:,1)+phi(:,:,2)+phi(:,:,3)).^erno)).^(1/erno);
max_error_bro     = max(error_bro(:));
max_error_brocor  = max(error_brocor(:));
mean_error_bro    = mean(error_bro(:));
mean_error_brocor = mean(error_brocor(:));
fprintf(1,['Maximum real space error: BRO=%.8g, BROCOR=%g (%2.3f%%)\n' ...
           'Mean real space error: BRO=%.8g, BROCOR=%g (%2.3f%%)\n'], ...
        max_error_bro,max_error_brocor,100*max_error_brocor/max_error_bro, ...
        mean_error_bro, mean_error_brocor,100*mean_error_brocor/ ...
        mean_error_bro);
ErrorMaps(error_bro,error_brocor,'Delta_Dot_');


