function [phase_map,intensity,phi1,phi2]=Beads(alpha,distance,energy,xshift,sign_of_phase,pow,alg)

if (nargin<1),alpha=12;end;
if (nargin<2),distance=1;end;
if (nargin<3),energy=30;end;
if (nargin<4),xshift=0;end;
if (nargin<5),sign_of_phase=-1;end;
if (nargin<6),pow=1;end;%m
if (nargin<7),alg=2;end;%keV

% Read data.
%distance  = 50; %m
slice     = 0; %first slice
%energy    = 10; %energy in keV
lambda    = EnergyConverter(energy);
pixelsize = 15e-7; %detector pixelsize in m
padding   = 1; %padded dim = padding*dim
padvalue  = 0;
iterations= 0;
%xshift    = 200;
% Data.
folder_name = sprintf('/home/moosmann/data/beads')';
file_name_prefix = sprintf('%s/beads_S%03u_E%02ukeV_R2048_P15e-7m_',folder_name,slice,energy);
file_name_postfix = 'usample-real.gpbin'; %unwrapped phase
[phase_map,phase_map,phase_map] ...
    = loadgpbin(sprintf('%s-%s',file_name_prefix,file_name_postfix));
[dimx,dimy] = size(phase_map);

phase_map = sign_of_phase*(phase_map - mean(phase_map(:)));
% Propagation.
[intensity,int_pad] = Propagation(phase_map,distance,lambda,pixelsize, ...
                                  padding);
% Reconstruction
padding = 1;
if alg==1,phi = reco(intensity,alpha,padding);
          phi = pixelsize^2/(2*pi*lambda*distance).*phi;
elseif alg==2,phi = Reco(intensity,alpha,lambda,distance,pixelsize,padding,padvalue,iterations);end;
phi1  = phi(:,:,1);
phi2  = phi(:,:,2);
domain(phase_map);
domain(phi1);
domain(phi2);
phi = phi1 + phi2;
domain(phi);
% Figures.
%figure('Name',sprintf('1D line cuts: alpha = %u, energy=%2.2u, d=%2.2u',alpha,energy,distance)),
%subplot(2,2,1),plot(phase_map(:,floor(dimy/2))),title('Exact phase map'),
%subplot(2,2,2),plot(phi (:,floor(dimy/2))),title('Bronnikov + Correction'),
%subplot(2,2,3),plot(phi1(:,floor(dimy/2))),title('Bronnikov'),
%subplot(2,2,4),plot(phi2(:,floor(dimy/2))),title('Correction');
x = 1:dimx;
figure('Name',sprintf(['1D line cuts: alpha = %u, energy=%2.2u, d=%u, black=exact_phase, ' ...
                    'blue=Bro, red=Bro+NLO'],alpha,energy,distance)),
plot(x,phase_map(:,xshift+floor(dimy/2)),'black',x,phi1(:,xshift+floor(dimy/2)),'blue', ...
     x,phi(:,xshift+floor(dimy/2)),'red');
saveas(gcf,sprintf('%s/beads_1d_linecuts.eps',folder_name),'psc2');
figure('Name',sprintf('2D phase maps: alpha = %u, energy=%2.2u, d=%u',alpha,energy,distance)),
subplot(2,2,1),imshow(phase_map,[]),colorbar,
subplot(2,2,2),imshow(phi, []),colorbar,
subplot(2,2,3),imshow(phi1,[]),colorbar,
subplot(2,2,4),imshow(phi2,[]),colorbar;
% Errors.
erb   = phase_map - phi1;
erbc  = phase_map - phi;
merb  = mean(abs(erb(:)).^pow).^(1/pow);
merbc = mean(abs(erbc(:)).^pow).^(1/pow);
fprintf(1,'Mean real space error: BRO=%.8g, BROCOR=%g (%2.3f%%)\n', ... 
        merb,merbc,100*merbc/merb);