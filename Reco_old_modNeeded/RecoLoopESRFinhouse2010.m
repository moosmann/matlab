function [phi1,phi2,intensity,absorption,phi] = RecoLoopESRFinhouse2010(images,scan_num,alpha,padding,include_absorption,iterations)

if nargin<1, images=0;end;
if nargin<2, scan_num = 3;end;
if nargin<3, alpha=4.2;end;
if nargin<4, padding=1;end;
if nargin<5, include_absorption=0;end;


write_edf = 0;
edf_crop = 769:1280;
x = 90:1909;
padvalue = 0;

% Read Data.
%scan_name = 'CT_meshesPApet_B_';
scan_name = 'CT_graphiteMinePApet_';
%scan_name = 'CT_carbonMeshesPApet_';
%scan_name = 'Stabheuschrecke_oben_';
scan_folder = sprintf('%s%u_',scan_name,scan_num);
file_name   = sprintf('%s%u_',scan_name,scan_num);
path_name   = sprintf('/mnt/tomoraid/tomo/ESRF_20100411_InhouseExperiment/CT_graphiteMinePApet/%s',scan_folder);
output_path = sprintf('/mnt/tomoraid/tomo/ESRF_20100411_InhouseExperiment/CT_graphiteMinePApet/%srecoBro_%u_/',scan_name,scan_num);
%path_name   = sprintf('/data/id19/inhouse/Lukas/%s',scan_folder);
%output_path = sprintf('/data/id19/inhouse1/Lukas/%srecoBro_%u_/',scan_name,scan_num);

[header,data] = pmedf_read(sprintf('%s/%s0000.edf',path_name,file_name));
flat1 = pmedfread(sprintf('%s/refHST0000.edf',path_name));
flat2 = pmedfread(sprintf('%s/refHST0500.edf',path_name));
flat3 = pmedfread(sprintf('%s/refHST1000.edf',path_name));
flat4 = pmedfread(sprintf('%s/refHST1500.edf',path_name));
dark  = pmedfread(sprintf('%s/darkend0000.edf',path_name));

% Extract Parameters.
energy    = str2num(pmedf_findInHeader(header,'energy')); % keV
lambda    = EnergyConverter(energy);% metre
pixelsize = str2num(pmedf_findInHeader(header,'optic_used')); % microns
pixelsize = pixelsize*1e-6; %metre
motor_pos = str2num(pmedf_findInHeader(header,'motor_pos')); 
distance  = motor_pos(10)*1e-3;% metre
fprintf(1,'energy=%g,  pixelsize=%g,  distance=%g\n',energy,pixelsize,distance);
if write_edf==1,
fprintf(1,'write edfs to: %s\n',output_path);
end;
% Loop.
for ii=images,
% Read data.
if ii>0,
data = pmedfread(sprintf('%s/%s%04u.edf',path_name,file_name,ii));
end;
% Choose mean flat field.
if ii<=250,                flat=flat1;
elseif (ii<=750)&(ii>250) ,flat=flat2;
elseif (ii<=1250)&(ii>750),flat=flat3;
elseif ii>1250,            flat=flat4;
end;
% Compute intensity.
intensity = (data(:,x) - dark(:,x))./(flat(:,x) - dark(:,x));
% Include absorption by division.
if (include_absorption&(scan_num>1)),
scan_folder = sprintf('%s%u_',scan_name,1);
file_name   = sprintf('%s%u_',scan_name,1);
path_name   = sprintf('/data/id19/inhouse1/Lukas/%s',scan_folder);
data  = pmedfread(sprintf('%s/%s%04u.edf',path_name,file_name,ii));
flat  = pmedfread(sprintf('%s/refHST0000.edf',path_name));
dark  = pmedfread(sprintf('%s/darkend0000.edf',path_name));
absorption = (data(:,x) - dark(:,x))./(flat(:,x) - dark(:,x));
intensity = intensity./absorption;
else, absorption = intensity;
end;
% Reconstruction.
phi = Reco(intensity,alpha,lambda,distance,pixelsize,padding,padvalue,iterations);
% Write retrieved phase to disc.
display(size(phi));
if write_edf,
edfwrite(sprintf('%s%s%04u_bro.edf',output_path,scan_name,ii),phi(:,edf_crop,1),'float32');
edfwrite(sprintf('%s%s%04u_cor.edf',output_path,scan_name,ii),phi(:,edf_crop,2),'float32');
end;
end;

% Rotate images.
phi1 = phi(:,:,1)';
phi2 = phi(:,:,2)';
intensity  = intensity';
absorption = absorption';
