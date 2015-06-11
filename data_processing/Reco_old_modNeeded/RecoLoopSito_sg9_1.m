function [phi,intensity] = RecoLoopSito_sg9_1(imfirst,imlast,crop,write_edf)

if nargin<3,crop = 0;end;
if nargin<4,write_edf = 0;end;

% Parameters.
%write_edf = 0;
alpha = 3.5; % SG9_1 3.61;%sg8;
energy    = 25; % keV
lambda    = EnergyConverter(energy);% metre
pixelsize = 6.6e-6; %metre
distance  = .4;% metre
padding = 1;
padvalue = 0;
iterations = 0;

% Directory and File names.
% INPUT folder:
folder=sprintf('/mnt/tomoraid3/tomo/TopoTomo-100420_InVivoTomography/Sitophilus_granarius/reco/sg9_1_5000fps/ffcorr');
file_prefix=sprintf('%s/ffcorr_radio_',folder);
% OUTPUT folder:
output_folder=sprintf('/mnt/tomoraid3/tomo/TopoTomo-100420_InVivoTomography/Sitophilus_granarius/reco/sg9_1_5000fps/phase_retrieval');

% Loop.
for ii=imfirst:imlast,

% Read intensity.
intensity = pmedfread(sprintf('%s%06g.edf',file_prefix,ii));
% Reconstruction.
if crop>0,intensity = intensity(crop:end,:);end;
[phi] = Reco(intensity,alpha,lambda,distance,pixelsize,padding,padvalue,iterations);
% Write retrieved phase to disc.
if write_edf,
    if (write_edf==1)&(ii==imfirst),fprintf(1,'write edfs to: %s\n',output_folder);end;
    edfwrite(sprintf('%s/Bro_%06u.edf',output_folder,ii),phi(:,:,1)','float32');
    edfwrite(sprintf('%s/Cor_%06u.edf',output_folder,ii),phi(:,:,2)','float32');
end;

end;

% Print parameters
fprintf(1,'energy = %gkeV, pixelsize = %gm, distance = %gm, resolution = %g x %g, regularisaion parameter = %g\n',energy,pixelsize,distance,size(intensity,1),size(intensity,2),alpha);

fprintf(1,'output format: %g x %g\n',size(intensity,1),size(intensity,2));
unix(sprintf('cp /home/moosmann/matlab/RecoLoopSito_sg9_1.m %s/',output_folder));
