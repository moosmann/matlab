function [phi1,phi2,intensity,absorption] = Stabheuschrecke(images,alpha,padding)

if nargin<1, images=0;end;
if nargin<2, alpha=4.2;end;
if nargin<3, padding=1;end;

write_edf = 0;
scan_num = 2;
x = 1:2048;

% Read Data.
scan_name   = 'Stabheuschrecke_oben_';
scan_folder = sprintf('%s%u_',scan_name,scan_num);
file_name   = sprintf('%s%u_',scan_name,scan_num);
path_name = sprintf('/data/id19/inhouse1/Lukas/%s',scan_folder);
ouput_path = path_name;

[header,data] = pmedf_read(sprintf('%s/%s0000.edf',path_name,file_name));
flat  = pmedfread(sprintf('%s/ref0000_0000.edf',path_name));
dark  = pmedfread(sprintf('%s/darkend0000.edf',path_name));

% Extract Parameters.
energy    = str2num(pmedf_findInHeader(header,'energy')); % keV
lambda    = EnergyConverter(energy);% metre
pixelsize = str2num(pmedf_findInHeader(header,'energy')); % microns
pixelsize = pixelsize*1e-6; %metre
motor_pos = str2num(pmedf_findInHeader(header,'motor_pos')); 
distance  = motor_pos(10)*1e-3;% metre
fprintf(1,'energy=%g,  pixelsize=%g,  distance=%g\n',energy,pixelsize,distance);


% Loop.
for ii=images,
if ii>0,
data = pmedfread(sprintf('%s/%s%04u.edf',path_name,file_name,ii));
end;
intensity = (data(x,:) - dark(x,:))./(flat(x,:) - dark(x,:));
if scan_num>1,
scan_folder = sprintf('%s%u_',scan_name,1);
file_name   = sprintf('%s%u_',scan_name,1);
path_name   = sprintf('/data/id19/inhouse1/Lukas/%s',scan_folder);
data  = pmedfread(sprintf('%s/%s%04u.edf',path_name,file_name,ii));
flat  = pmedfread(sprintf('%s/ref0000_0000.edf',path_name));
dark  = pmedfread(sprintf('%s/darkend0000.edf',path_name));
absorption = (data(x,:) - dark(x,:))./(flat(x,:) - dark(x,:));
intensity = intensity./absorption;
end;

phi = Reco(intensity,alpha,lambda,distance,pixelsize,padding);
if write_edf,
edfwrite(sprintf('%s%s%04u_bro.edf',output_path,folder_name,ii),phi(:,:,1),'float32');
edfwrite(sprintf('%s%s%04u_cor.edf',output_path,folder_name,ii),phi(:,:,2),'float32');
end;
end;

% Rotate images.
phi1 = phi(:,:,1)';
phi2 = phi(:,:,2)';
intensity  = intensity';
absorption = absorption';
