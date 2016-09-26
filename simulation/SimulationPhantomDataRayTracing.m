function [phase_map,intensity,phi,intensity2,phi2] = Simulation(energy,distance,binning,hsize,sigma,print,shifted_plots,binned_propagation,raytracing)
% Simulation of phantom data via raytracing (SRCLsim), propagation,
% binning and phase retrieval.
% Input: enery in keV, distance in m, raytracing yes (1,default)/no (0)

    if nargin<1, energy   = 30;end; % kev!!
    if nargin<2, distance = 0.5;end; % metres!!
    if nargin<3, hsize = 1;end;
    if nargin<4, sigma = 0;end;
    if nargin<5, print = 0;end;
    if nargin<6, shifted_plots = 1;end;
    if nargin<7, binned_propagation = 1;end;
    if nargin<8, raytracing = 0;end;
    
% 'PMMA' at 10keV: 2.67147e-06, at 30keV: 2.96111e-07
% 'protein' at 10keV: 2.99405e-06, protein at 30keV: 3.31785e-07

% PARAMETERS.
compound      = 'PMMA';
slice         = 0;
lambda        = EnergyConverter(energy);
ddelta        = 0.05;
object_folder = '/home/moosmann/data/raytracing';
detector_size = 0.002048;
% Before binning: resolution = bin_fac^2 x (2048 x 2048), detector size =
% 0.002048m -> pixelsize = detector_size/resolution = 0.001m/bin_fac
alpha      = 12;
padding    = 1;
padvalue   = 1;
iterations = 0;
pow        = 1;
print_intensity = 0;

% Load SRCLsim links into bash and compute line integrals along the ray path
% for every pixel using SRCLsim.
eval(sprintf('cd %s',object_folder));
if raytracing,
tic;
unix('source ~/.bashrc');
unix('. ~/bin/SRCLsim/bin/setup_SRCLsim_bash.sh');
unix(sprintf('SRCLsim-SamplePhaseShift  %s/setup.srcl',object_folder));
traytracing = toc;
else, traytracing = 0;
end;

% Read object into Matlab file.
tic;
object_file = sprintf('%s/half_bead_%04u.edf',object_folder,slice);
phase_map = pmedfread(object_file);
treading = toc;
[resx_unbinned,resy_unbinned] = size(phase_map);
resx_binned = resx_unbinned/binning;
resy_binned = resx_binned;
pixelsize_unbinned = detector_size/resx_unbinned;
pixelsize_binned   = detector_size/resx_binned;
binx               = resx_unbinned/resx_binned;
biny               = resy_unbinned/resy_binned;

% Compute Refractive Index.
tic;
ref_ind_file = sprintf('%s/ref_ind.txt',object_folder);
unix(sprintf('pmasf -v -e%u -C %s | grep delta | sed ''s/delta = 1-n = (//'' | sed ''s/,/ /'' | sed ''s/)//'' | sed ''s/\\W//'' > %s',1000*energy,compound,ref_ind_file));
[delta,beta] = textread(ref_ind_file);
tpmasf = toc;
fprintf('Compound: %s, energy: %gkeV, refractive index: [delta,beta]=[%g,%g]\n',compound,energy,delta,beta);
% Create Phantom.
phase_map = 2*pi/lambda*delta*(phase_map + (1+ddelta)*flipud(phase_map));
phase_map = phase_map - mean(phase_map(:));
% Blurring.
if sigma>0,
    phase_map = imfilter(phase_map,fspecial('gaussian',[hsize hsize],sigma));
    phase_map = phase_map - mean(phase_map(:));
end;
%%%%%%%%%%%%%%%%%% Propagation of unbinned phase map. %%%%%%%%%%%%%%%%%%%%%%%
tic;
fprintf('\nPROPAGATION OF UNBINNED PHASE MAP\n');
tprop = toc;
[intensity] = Propagation(phase_map,distance,lambda,pixelsize_unbinned,padding);
%%%%%%%% CHANGE OBJECT_FOLDER %%%%%%%%%%%
object_folder_main = object_folder;
object_folder = sprintf('%s/%s_E%ukev_d%ucm_%ux%u',object_folder,compound,energy,100*distance,resx_unbinned,resy_unbinned);
unix(sprintf('mkdir %s',object_folder));
% Binning Of Intensity.
tic;
%figure('Name','Intensity unbinned. Propagation unbinned.'),imshow(intensity,[]),colorbar;
[intensity] = Binning(intensity,[resx_binned resy_binned]);
edfwrite(sprintf('%s/intensity_BinningAfterPropagation.edf',object_folder),intensity,'float32');
%figure('Name','Intensity binned. Binning after propagation'),imshow(intensity,[]),colorbar;
% Binning of phase_map.
[phase_map] = Binning(phase_map,[resx_binned resy_binned]);
tbinning = toc;
% Phase retrieval.
tic;
[phi]  = Reco(intensity,alpha,lambda,distance,pixelsize_binned,padding,padvalue,iterations);
treco = toc;
% Domain.
Domain(phi(:,:,1),'retrieved phase');
Domain(phi(:,:,2),'retrieved correction');
% Real space error measures.
[merb,merbc,sterb,sterbc,erb,erbc] = Errors(phase_map,phi,pow);
error_file = sprintf('%s/error.txt',object_folder);
fid = fopen(error_file,'wt');
fprintf(fid,['Real space errors. Binning after Propagation: MeanBro=%g, ' ...
           'MeanBro=%g (%2.3f%%), StaBro=%g, StaBroCor=%g\n'],merb,merbc,100*merbc/merb,sterb,sterbc);
%%%%%%%%%%%%%%%%% Propagation after Binning. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if binned_propagation,
fprintf('\nPROPAGATION OF BINNED PHASE MAP\n')
[intensity2] = Propagation(phase_map,distance,lambda,pixelsize_binned,padding);
edfwrite(sprintf('%s/intensity_BinningBeforePropagation.edf',object_folder),intensity2,'float32');
% Phase retrieval.
tic;
[phi2] = Reco(intensity2,alpha,lambda,distance,pixelsize_binned,padding,padvalue,iterations);    
treco = treco + toc;
% Print computational time, error measures and figures.
%fprintf(1,['Duration: raytracing: %gs, reading: %gs, pmasf: %gs, propagation: %gs, binning: %gs, phase retrieval: %gs\n'],traytracing,treading,tpmasf,tprop,tbinning,treco);
% Domain.
Domain(phi2(:,:,1),'retrieved phase');
Domain(phi2(:,:,2),'retrieved correction');
% Real space error measures.
[merb,merbc,sterb,sterbc,erb,erbc] = Errors(phase_map,phi2,pow);
fprintf(fid,['Real space errors. Binning after Propagation: MeanBro=%g, ' ...
           'MeanBro=%g (%2.3f%%), StaBro=%g, StaBroCor=%g\n'],merb,merbc,100*merbc/merb,sterb,sterbc);
end;
fclose(fid);

% Write files to harddrive.
edfwrite(sprintf('%s/phase_rec_bro.edf',object_folder),phi(:,:,1),'float32');
edfwrite(sprintf('%s/phase_rec_cor.edf',object_folder),phi(:,:,2),'float32');

%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = ceil(resx_binned/4):floor(resx_binned*3/4);
y = floor(resy_binned/2);
if print,
% Retrieved phase plots.
LinePlots(phase_map,phi,sprintf('Retrieved phase. Binning after propagation (E=%gkeV,d=%g,sigma=%g):',energy,distance),y,x),
saveas(gcf,sprintf('%s/LineCut_y1d2_RetrievedPhase_BinAfterProp.eps',object_folder),'psc2');
saveas(gcf,sprintf('%s/LineCut_y1d2_RetrievedPhase_BinAfterProp.tiff',object_folder),'tiffn');
if binned_propagation,
LinePlots(phase_map,phi2,sprintf('Retrieved phase. Binning before propagation (E=%gkeV,d=%g,sigma=%g):',energy,distance),y,x);end;
%%%%%%%% Y SHIFTED LINE CUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
if shifted_plots,
y = floor(resy_binned*(19/32));
% Retrieved phase plots.
LinePlots(phase_map,phi,sprintf('Retrieved phase. Binning after propagation (E=%gkeV,d=%g):',energy,distance),y,x),
saveas(gcf,sprintf('%s/LineCut_y7d16_RetrievedPhase_BinAfterProp.eps',object_folder),'psc2');
saveas(gcf,sprintf('%s/LineCut_y7d16_RetrievedPhase_BinAfterProp.tiff',object_folder),'tiffn');
if binned_propagation,
LinePlots(phase_map,phi2,sprintf('Retrieved phase. Binning before propagation (E=%gkeV,d=%g):',energy,distance),y,x);end;
end;
end;

unix(sprintf('cp %s/ref_ind.txt %s/',object_folder_main,object_folder));
unix(sprintf('cp %s/setup.srcl %s/',object_folder_main,object_folder));
unix(sprintf('cp %s/sample.srcl %s/',object_folder_main,object_folder));

fprintf('\n');