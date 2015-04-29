function PaperIntPhase(phex,int,phi,distance);

if nargin<4,distance=1;end;

    
resolution = size(phex,1);
folder = '/home/moosmann/data/SiemensStar/paper_figures';
file_name_prefix = sprintf('%s/Map2D_E30keV_z%03gcm_res%04g',folder,100*distance,resolution);

fip=figure('color','white','Name','intensity and exact phase map');
imshow([phex,int-1],[]);colormap(flipud(colormap('gray')));colorbar;
set(gca,'XTick',[],'YTick',[]);
flonlo=figure('color','white','Name','reconstructed phase + correction map');
imshow([phi(:,:,1),phi(:,:,2)],[]);colormap(flipud(colormap('gray')));colorbar;
set(gca,'XTick',[],'YTick',[]);
fer=figure('color','white','Name','error map');
imshow([abs(phi(:,:,1)-phex),abs(phi(:,:,1)+phi(:,:,2)-phex)],[]);colormap(flipud(colormap('gray')));colorbar;
set(gca,'XTick',[],'YTick',[]);

saveas(fip,sprintf('%s_intensity_phase.eps',file_name_prefix,'epsc'));
saveas(flonlo,sprintf('%s_phaseNL_phaseNLO.eps',file_name_prefix,'epsc'));
saveas(fer,sprintf('%s_errorNL_errorNLO.eps',file_name_prefix,'epsc'));
