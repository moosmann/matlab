function [phi,intensity,flat] = RecoLoopSito_sg9_1(imfirst,imlast,flat_end,x_first,write_edf)

if nargin<3,flat_end = 200;end;
if nargin<4,x_first = 1;end;
if nargin<5,write_edf = 0;end;

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
%folder=sprintf('/mnt/tomoraid3/tomo/TopoTomo-100420_InVivoTomography/Sitophilus_granarius/hs-tomo/sg8_5000fps_1');
% INPUT folder:
folder=sprintf('/mnt/tomoraid3/tomo/TopoTomo-100420_InVivoTomography/Sitophilus_granarius/hs-tomo/sg9_1_5000fps');
file_prefix=sprintf('%s/radio_C001H001S0001',folder);
% OUTPUT folder:
output_folder=sprintf(['/mnt/tomoraid3/tomo/TopoTomo-100420_InVivoTomography/Sitophilus_granarius/reco/sg9_1_5000fps/phase_retrieval']);
% Read flat fields.
flat=double(imread(sprintf('%s%06g.tif',file_prefix,1)));
for ii=2:flat_end, flat=flat+double(imread(sprintf('%s%06g.tif',file_prefix,ii)));end;
flat = flat/flat_end;
% Print parameters
fprintf(1,'energy = %gkeV, pixelsize = %gm, distance = %gm, resolution = %g x %g\n',  ...
 energy,pixelsize,distance,size(flat,1),size(flat,2));

% Loop.
for ii=imfirst:imlast,

% Read intensity.
intensity = double(imread(sprintf('%s%06g.tif',file_prefix,ii)))./flat;
% Reconstruction.
intensity = intensity(x_first:end,:);
[phi] = Reco(intensity,alpha,lambda,distance,pixelsize,padding,padvalue,iterations);
% Write retrieved phase to disc.
if write_edf,
    if (write_edf==1)&(ii==imfirst),fprintf(1,'write edfs to: %s\n',output_folder);end;
    edfwrite(sprintf('%s/Bro%06u.edf',output_folder,ii),phi(:,:,1)','float32');
    edfwrite(sprintf('%s/Cor%06u.edf',output_folder,ii),phi(:,:,2)','float32');
end;

end;

fprintf(1,'output format: %g x %g\n',size(intensity,1),size(intensity,2));
%unix(sprintf('cd %s',output_folder)); %not working in function body
unix(sprintf('cp /home/moosmann/matlab/RecoLoopSito_sg9_1 %s/',output_folder));
