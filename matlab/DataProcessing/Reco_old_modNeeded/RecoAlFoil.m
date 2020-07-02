function [im,data,ref,dark,darkend,LO,NLO]=RecoAlFoil(doReco,alpha,kend)
% Phase retrieval to LO and NLO for data from Experiment at
% ID22@ESRF

if nargin<2, alpha = 4;end;
if nargin<3, kend = 4; end;

% Reco.m arguments.
% ID22@ESRF
energy     = 17.5;
lambda     = EnergyConverter(energy);
distance   = 1.1590;
pixelsize  = 0.0606e-6;
padvalue   = 0;
iterations = 0;
compute_correction = 1;
padding    = 1 ;

% Folder to data.
ParentFolder = '/mnt/tomoraid3/user/moosmann/Al_foil_B_1_/';
dark         = edfread([ParentFolder 'dark.edf']);
darkend      = edfread([ParentFolder 'darkend0000.edf']);

% Loop over data.
for kk=1:kend,
    proj = [0 500 1000 1500];
            dat          = edfread([ParentFolder 'Al_foil_B_1_' num2str(proj(kk),'%04u') '.edf']);
            refHST       = edfread([ParentFolder 'refHST' num2str(proj(kk),'%04u') '.edf']);
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
