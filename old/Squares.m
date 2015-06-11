function [phase_map,intensity,phi1,phi2]=Squares(distance)


alpha=12;
alg=2;
pow=1;

% Read data.
slice     = 0; %first slice
energy    = 30; %energy in keV
lambda    = EnergyConverter(energy);
pixelsize = 15e-7; %detector pixelsize in m
padding   = 1; %padded dim = padding*dim
% Data.
cd ~/data/squares;
file_name_prefix = sprintf('squares_S%03u_E%02ukeV_R2048_P15e-7m_',slice,energy);
file_name_postfix = 'usample-real.gpbin'; %unwrapped phase
[phase_map,phase_map,phase_map] ...
    = loadgpbin(sprintf('%s-%s',file_name_prefix,file_name_postfix));
[dimx,dimy] = size(phase_map);
phase_map = -(phase_map - mean(phase_map(:)));
% Propagation.
[intensity,int_pad] = Propagation(phase_map,distance,lambda,pixelsize, ...
                                  padding);
% Reconstruction
padding = 1;
if alg==1,phi = reco(intensity,alpha,padding);
          phi = pixelsize^2/(2*pi*lambda*distance).*phi;
elseif alg==2,phi = Reco(intensity,alpha,lambda,distance,pixelsize,padding);end;
phi1  = phi(:,:,1);
phi2  = phi(:,:,2);
domain(phase_map);
domain(phi1);
domain(phi2);
phi = phi1 + phi2;
domain(phi);
% Figures.
figure('Name',sprintf('alpha = %u',alpha)),
subplot(2,2,1),plot(phase_map(:,floor(dimy/2))),title('Exact phase map'),
subplot(2,2,2),plot(phi (:,floor(dimy/2))),title('Bronnikov + Correction'),
subplot(2,2,3),plot(phi1(:,floor(dimy/2))),title('Bronnikov'),
subplot(2,2,4),plot(phi2(:,floor(dimy/2))),title('Correction');
x = 1:dimx;
figure('Name',sprintf('alpha = %u',alpha)),
plot(x,phase_map(:,floor(dimy/2)),'black',x,phi1(:,floor(dimy/2)),'blue',x,phi(:,floor(dimy/2)),'red');
figure('Name',sprintf('alpha = %u',alpha)),
subplot(2,2,1),imshow(phase_map,[]),colorbar,
subplot(2,2,2),imshow(phi, []),colorbar,
subplot(2,2,3),imshow(phi1,[]),colorbar,
subplot(2,2,4),imshow(phi2,[]),colorbar;
% Errors.
erb   = phase_map - phi1;
erbc  = phase_map - phi;
merb  = mean(abs(erb(:)).^pow).^-pow;
merbc = mean(abs(erbc(:)).^pow).^-pow;
fprintf(1,'Mean real space error: BRO=%.8g, BROCOR=%g (%2.3f%%)\n', ... 
        merb,merbc,100*merbc/merb);
