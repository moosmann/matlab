function [phi,intensity] = RecoLoopBug(imfirst,imlast,x_first,write_edf)

if nargin<3,x_first = 0;end;
if nargin<4,write_edf = 0;end;

% Parameters.
%write_edf = 0;
alpha = 4.7; % from Lcurve
energy    = 20; % keV pink beam
lambda    = EnergyConverter(energy);% metre
pixelsize = 0.645e-6; %metre
distance  = .087;% metre
padding = 1;
padvalue = 0;
iterations = 0;

% Directory and File names.
% INPUT FOLDER: 
folder=sprintf('/mnt/tomoraid3/user/moosmann/TopoTomo-100420-BiologicalSamples/analysis_pv/Cosmo-complete/radiograms');
file_prefix=sprintf('%s/radio_',folder);
% OUTPUT FOLDER:
output_folder=sprintf(['/mnt/tomoraid3/user/moosmann/TopoTomo-100420-BiologicalSamples/analysis_pv/Cosmo-complete/phase_retrieval']);

% Loop.
for ii=imfirst:imlast,

% Read intensity.
intensity = double(imread(sprintf('%s%04g.tif',file_prefix,1)));
% Reconstruction.
if x_first>0,intensity = intensity(x_first:end,:);end;
[phi] = Reco(intensity,alpha,lambda,distance,pixelsize,padding,padvalue,iterations);
% Write retrieved phase to disc.
if write_edf,
    if (write_edf==1)&(ii==imfirst),fprintf(1,'write edfs to: %s\n',output_folder);end;
    edfwrite(sprintf('%s/Bro_%04u.edf',output_folder,ii),phi(:,:,1)','float32');
    edfwrite(sprintf('%s/Cor_%04u.edf',output_folder,ii),phi(:,:,2)','float32');
end;

end;

% Print parameters
fprintf(1,'energy = %gkeV, pixelsize = %gm, distance = %gm, resolution = %g x %g\n',  ...
 energy,pixelsize,distance,size(intensity,1),size(intensity,2));
fprintf(1,'output format: %g x %g\n',size(intensity,1),size(intensity,2));
%unix(sprintf('cd %s',output_folder)); %not working in function body
unix(sprintf('cp /home/moosmann/matlab/RecoLoopBug.m %s/',output_folder));
