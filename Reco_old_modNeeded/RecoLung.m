function  [data,phaseCl,LO,NLO]=RecoLung(doReco,alpha)
% Phase retrieval to LO and NLO for data from Inhouse Experiment at
% ID19@ESRF in April

if nargin<2, alpha = 0;end;

% Reco.m arguments.
% ESRF ID22 experimental setup parameters.
energy     = 17.3;
lambda     = EnergyConverter(energy);
distance   = [0.030061   0.031025   0.034490   0.042886];
pixelsize  = 5.8539e-08;
if alpha == 0,
    alpha = [4.9 4.83 4.76 4.79];
end;
padvalue   = 0;
iterations = 0;
compute_correction = 1;
padding    = 1 ;

% Folder to data.
ParentFolder = '/mnt/tomoraid3/user/moosmann/Lung_sample__Suhonnen/';
phaseCl        = edfread([ParentFolder 'lung_sample2_highres_0000.edf']);

% Loop over data.
for kk=1:4,
            dat          = edfread([ParentFolder 'lung_sample2_highres1_' num2str(kk) '_flt0000.edf']);
            data(:,:,kk) = dat;
            % Reconstruct phase
            if doReco,
                if ndims(alpha) > 1,
                    phi = Reco(dat,alpha(kk),lambda,distance(kk),pixelsize, ...
                           padding,padvalue,iterations, ...
                               compute_correction);
                else,
                    phi = Reco(dat,alpha,lambda,distance(kk), ...
                               pixelsize,padding,padvalue,iterations,compute_correction);
                end;
                
            LO(:,:,kk)   = phi(:,:,1);
            NLO(:,:,kk)  = phi(:,:,2);
            edfwrite([ParentFolder 'phaseLO_' num2str(kk) '.edf'],phi(:,:,1),'float32');
            edfwrite([ParentFolder 'phaseNLO_' num2str(kk) '.edf'],NLO,'float32');
            LinePlots(phaseCl,phi,'',1000); 
            end;
end;
