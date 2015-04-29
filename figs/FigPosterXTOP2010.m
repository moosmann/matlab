%function PosterPics(hsize,sigma,dist)

energy      = 30; %keV
%dist        = 0.3; %m
resolution  = 2048;
spokes      = 256;
binfac      = 1;
hsize       = 8;
sigma       = 1.5;
dist        = 0.3;
delta       = 1e-7;
linecutat   = floor(0.4395*resolution);

[pha,int,phi]=SimulationSiemensStar(energy,dist,resolution,binfac,hsize,sigma,spokes,delta);

picfolder = '/home/moosmann/poster/';
prefix    = sprintf(['%sE%2ikeV_z%2icm_'],picfolder,energy,100*dist);
fprintf(1,['Prefix: ' prefix '\n']);
fprintf(1,['Line cut at ' num2str(linecutat) '\n']);

pha = -pha;
lo = -phi(:,:,1);
nlo = -phi(:,:,2);

erlo = -abs(lo-pha);
ernlo = -abs(lo+nlo-pha);


edfwrite([ prefix 'ExactPhase.edf'],pha,'float32');
edfwrite([ prefix 'PhaseLO.edf'],lo,'float32');
edfwrite([ prefix 'PhaseNLO.edf'],nlo,'float32');
edfwrite([ prefix 'ErrorLO.edf'],erlo,'float32');
edfwrite([ prefix 'ErrorNLO.edf'],ernlo,'float32');





