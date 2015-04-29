function [phase_map,intensity,phi,partab,intensity2,phi2] = SimulationSiemensStar(energy,distance,resolution,binning_factor,hsize,sigma,spokes,delta,noise,print_oversampling)
% Simulation of Siemens star pattern, propagation, binning and phase
% retrieval.  Input: enery in keV, distance in m, raytracing yes
% (1,default)/no (0)

if nargin<1, energy   = 30;end; % kev!!
if nargin<2, distance = 0.5;end; % metres!!
if nargin<3, resolution = 512;end;
if nargin<4, binning_factor = 4;end;
if nargin<5, hsize = 5;end;
if nargin<6, sigma = 0;end;
if nargin<7, spokes = 0;end;
if nargin<8, delta = 1e-7;end;
if nargin<9, noise = 0;end;
if nargin<10, print_oversampling = 0;end;

% PARAMETERS.
lambda        = EnergyConverter(energy);
object_folder = '/home/moosmann/data/SiemensStar';
error_file    = sprintf('%s/error.txt',object_folder);
par_file      = sprintf('%s/parameters.txt',object_folder);
%resx_binned   = 2048;
%resy_binned   = 2048;
detector_size = 0.002048;
thickness     = detector_size/8;
% Before binning: resolution = bin_fac^2 x (2048 x 2048), detector size =
% 0.002048m -> pixelsize = detector_size/resolution = 0.001m/bin_fac
alpha         = 12;
padding       = 1;
padvalue      = 1;
iterations    = 0;
pow           = 1;
print_time    = 0;
unbinned_propagation = 0;

% Siemens star pattern.
star               = SiemensStar(resolution,spokes);
if sigma>0,star    = imfilter(star,fspecial('gaussian',[hsize hsize],sigma));end;
phase_map          = 2*pi/lambda*delta*thickness*star;
phase_map          = phase_map - mean(phase_map(:));
[resx_unbinned,resy_unbinned] = size(phase_map);
resx_binned        = resx_unbinned/binning_factor;
resy_binned        = resy_unbinned/binning_factor;
pixelsize_unbinned = detector_size/resx_unbinned;
pixelsize_binned   = detector_size/resx_binned;
binx               = resx_unbinned/resx_binned;
biny               = resy_unbinned/resy_binned;

%%%%%%%%%%%%%%%%%% Propagation of unbinned phase map (before Binning).
fprintf('PROPAGATION OF UNBINNED PHASE MAP\n');
tic;
[intensity] = Propagation(phase_map,distance,lambda,pixelsize_unbinned,padding);
tprop       = toc;
% Add noise.
if noise>0,
    intensity   = imnoise(intensity,'poisson');
    Domain(intensity,'intensity with noise');
end;
% Binning Of Intensity.
tic;
[intensity] = Binning(intensity,[resx_binned resy_binned]);
tbinning    = toc;
% Phase retrieval.
tic;
[phi] = Reco(intensity,alpha,lambda,distance,pixelsize_binned,padding,padvalue,iterations);
treco = toc;
% Print computational time, error measures and figures.
if print_time,fprintf(1,['propagation: %gs, binning: %gs, phase retrieval: %gs\n'],tprop,tbinning,treco);end;
% Print domain of phase map.
[phex_min,phex_max]=Domain(phase_map,'exact phase');
% Binning of phase_map.
[phase_map] = Binning(phase_map,[resx_binned resy_binned]);
% Domain.
[phi1_min,phi1_max]=Domain(phi(:,:,1),'retrieved phase');
[phi2_min,phi2_max]=Domain(phi(:,:,1)+phi(:,:,2),'corrected phase');
% Real space error measures.
[merb,merbc,sterb,sterbc,erb,erbc] = Errors(phase_map,phi,pow);
fid    = fopen(error_file,'wt');
fprintf(fid,['Real space errors. Binning after Propagation: MeanBro=%g, ' ...
    'MeanBro=%g (%2.3f%%), StaBro=%g, StaBroCor=%g\n'],merb,merbc,100*merbc/merb,sterb,sterbc);
par    = fopen(par_file,'a');
partab = sprintf( ['%3g%5g%9.4g%9.3g%5g' ...
    '%7.0u%7.0u%8.0g%8.0g%5.0u' ...
    '%5g%5g%5g%8.4f%8.4f' ...
    '%8.4f%8.4f%8.4f%8.4f%12.3e' ...
    '%12.3e%6.3g%%%12.3e%12.3e\n'], ...
    energy,distance,detector_size,delta,spokes, ...
    resx_binned,resx_unbinned,pixelsize_binned,pixelsize_unbinned,alpha, ...
    sigma,hsize,binning_factor, phex_min,phex_max, ...
    phex_max-phex_min,phi1_min,phi1_max,phi1_max-phi1_min,merb, ...
    merbc,100*merbc/merb,sterb,sterbc);
fprintf(par,partab);

%%%%%%%%%%%%%%%%% Propagation after Binning.
if unbinned_propagation(1)
    fprintf('\nPROPAGA TION OF BINNED PHASE MAP\n')
    tic;
    [intensity2] = Propagation(phase_map,distance,lambda,pixelsize_binned,padding);
    tprop        = toc;
    % Add noise.
    if noise>0,
        intensity   = imnoise(intensity,'poisson');end;
    % Phase retrieval.
    tic;
    [phi2] = Reco(intensity2,alpha,lambda,distance,pixelsize_binned,padding,padvalue,iterations);
    treco = toc;
    % Print computational time, error measures and figures.
    if print_time,fprintf(1,['propagation: %gs, binning: %gs, phase retrieval: %gs\n'],tprop,tbinning,treco);end;
    % Domain.
    [phi1_min,phi1_max]=Domain(phi(:,:,1),'retrieved phase');
    [~,phi2_max]=Domain(phi(:,:,1)+phi(:,:,2),'corrected phase');
    % Real space error measures.
    [merb,merbc,sterb,sterbc,erb2,erbc2] = Errors(phase_map,phi2,pow);
    fprintf(fid,['Real space errors. Binning before Propagation: MeanBro=%g, ' ...
        'MeanBro=%g (%2.3f%%), StaBro=%g, StaBroCor=%g\n\n'],merb,merbc,100*merbc/merb,sterb,sterbc);
    %edfwrite(sprintf('%s/intensity_BinningBeforePropagation.edf',object_folder),intensity2,'float32');
    fprintf(par, ['%3g%5g%9.4g%9.3g%5g' ...
        '%7.0u%7.0u%8.0g%8.0g%5.0u' ...
        '%5g%5g%5g%8.4f%8.4f' ...
        '%8.4f%8.4f%8.4f%8.4f%12.3e' ...
        '%12.3e%6.3g%%%12.3e%12.3e\n'], ...
        energy,distance,detector_size,delta,spokes, ...
        resx_binned,resx_unbinned,pixelsize_binned,pixelsize_unbinned,alpha, ...
        sigma,hsize,binning_factor, phex_min,phex_max, ...
        phex_max-phex_min,phi1_min,phi1_max,phi1_max-phi1_min,merb, ...
        merbc,100*merbc/merb,sterb,sterbc);
    
end;

% Close text files.
fclose(fid);
fclose(par);

%%%%%%%%%%%%%%%%
% Save images as edfs to harddrive.
edfwrite(sprintf('%s/intensity_BinningAfterPropagation.edf',object_folder),intensity,'float32');
edfwrite(sprintf('%s/phase_rec_bro.edf',object_folder),phi(:,:,1),'float32');
edfwrite(sprintf('%s/phase_rec_cor.edf',object_folder),phi(:,:,2),'float32');

% Intensity plots.
if print_oversampling,
    ishowcb(intensity,  'Intensity binned. Binning after propagation');
    ishowcb(phi(:,:,1), 'Phase map. Binning after propagation');
    ishowcb(erb,  'Error map. Binning before propagation: Bronnikov');
    ishowcb(erbc, ['Error map. Binning before propagation: Bronnikov + Correction']);
    if binned_propagation,
        ishowcb(intensity2, 'Intensity binned. Binning before propagation');
        ishowcb(phi2(:,:,1),'Phase map. Binning before propagation');
        ishowcb(erb2, 'Error map. Binning after propagation: Bronnikov');
        ishowcb(erbc2,['Error map. Binning after propagation: Bronnikov + Correction']);
    end;
end;

ishowcb([phi(:,:,1),phi(:,:,2)],sprintf(['Reconstructed phase maps (E=%gkeV,d=%g): Left: ' ...
    'Bronnikov, Right: Correction'],energy,distance));
ishowcb([abs(phase_map-phi(:,:,1)),abs(phase_map-phi(:,:,1)-phi(:,:,2))], ...
    sprintf('Absolute phase map errors (E=%gkeV,d=%g): Left: Bronnikov-Exact, Right: Bronnikov+Correction-Exact',energy,distance));

LinePlots(phase_map,phi,sprintf('Retrieved phase maps (E=%gkeV,d=%g):',energy,distance),floor(0.4395*resolution),1:2000);

fprintf('\n');
