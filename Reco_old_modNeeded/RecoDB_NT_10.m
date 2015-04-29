function  [data,LO,NLO]=RecoDB_NT_10(doReco,alpha)
% Phase retrieval to LO and NLO for data from Experiment at ESRF:
% aluminum block

if nargin<2, alpha = 0;end;

% Reco.m arguments.
% ESRF ID22 experimental setup parameters.
energy     = 20;
lambda     = EnergyConverter(energy);
distance   = 0.008;
pixelsize  = 0.7e-06;
padvalue   = 0;
iterations = 0;
compute_correction = 1;
padding    = 1 ;

% Folder to data.
ParentFolder = '/mnt/tomoraid3/user/moosmann/DB_NT_10_ffCorr/';

data          = edfread([ParentFolder 'proj0004.edf']);
% Reconstruct phase
if doReco,
    phi = Reco(data,alpha,lambda,distance, ...
                               pixelsize,padding,padvalue,iterations, ...
               compute_correction);
    LO = phi(:,:,1);
    NLO = phi(:,:,2);
 end;
