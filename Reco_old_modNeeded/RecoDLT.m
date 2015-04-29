function [im,data,ref,dark,LO,NLO]=RecoDLT(doReco,alpha,kend)
% Phase retrieval to LO and NLO for data from Experiment at
% ID19@ESRF (/data/id19/inhouse/sep00/fleece/6_1_)

% resolution: 1024 x 1024
% projections: 900
% magnification (optics): 6

if nargin<2, alpha = 4;end;
if nargin<3, kend = 4; end;

% Reco.m arguments.
% ID22@ESRF
energy     = 25; % lambda = 4.9594e-11 m
lambda     = EnergyConverter(energy);
distance   = 0.973;
pixelsize  = 6.7e-6;
padvalue   = 0;
iterations = 0;
compute_correction = 0;
padding    = 1 ;
ymin       = 151;
ymax       = 918;

% Folder to data.
ParentFolder = '/mnt/tomoraid3/user/moosmann/DLT/6_3_/';
dark         = edfread([ParentFolder 'dark.edf']);
dark         = dark(:,ymin:ymax);
%darkend      = edfread([ParentFolder 'darkend0000.edf']);

% Loop over data.
for kk=1:kend,
    proj = [0 300 600 900];
            dat          = edfread([ParentFolder '6_3_' num2str(proj(kk),'%04u') '.edf']);
            dat          = dat(:,ymin:ymax);
            refHST       = edfread([ParentFolder 'refHST' num2str(proj(kk),'%04u') '.edf']);
            refHST       = refHST(:,ymin:ymax);
            data(:,:,kk) = dat;
            ref(:,:,kk)  = refHST;
            im(:,:,kk)   = (dat-dark)./(refHST-dark);
            % Reconstruct phase
            if doReco,
                    phi = Reco(im(:,:,kk),alpha,lambda,distance(kk), ...
                               pixelsize,padding,padvalue,iterations,compute_correction);
            LO(:,:,kk)   = phi(:,:,1);
            NLO(:,:,kk)  = phi(:,:,2);
            %            edfwrite([ParentFolder 'phaseLO_' num2str(kk) '.edf'],phi(:,:,1),'float32');
            %            edfwrite([ParentFolder 'phaseNLO_' num2str(kk) '.edf'],NLO,'float32');
            end;
end;
