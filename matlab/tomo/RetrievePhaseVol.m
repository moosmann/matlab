function delta = RetrievePhaseVol(vol,Method,RegPar,RamLak,padAll,EnergyDistancePixelsize)

%% Defaults
if nargin < 2
    Method = 'quasi';
end
if nargin < 3
    RegPar = 2;
end
if nargin < 4
    RamLak = 0;
end
if nargin < 5
    padAll = 0;
end
if nargin < 6
    %EnergyDistancePixelsize = [20,0.945,0.75e-6];
    EnergyDistancePixelsize = [30, 0.704, 1.1e-6];
end

%% MAIN
tic
% padarray along 3rd dimension to avoid jumps of data at boundary which
% introduces ghost arifacts: first and last slices are superimposed
if padAll
    hpad = @(a) padarray(a,[size(a)/2],'symmetric','both');
    % crop array back to original dimension
    hcrop = @(a) a(size(a,1)/4+(1:size(a,1)/2), size(a,2)/4+(1:size(a,2)/2), size(a,3)/4+(1:size(a,3)/2));
else
    hpad = @(a) padarray(a,[0,0,size(a,3)/2],'symmetric','both');
    % crop array back to original dimension
    hcrop = @(a) a(:,:,size(a,3)/4+(1:size(a,3)/2));
end

% backprojection filter: 'Ram-Lak'
% if RamLak
%     RamLak = @(VolSize) repmat(sqrt(repmat(FrequencyVector(VolSize(2)),[VolSize(1) 1]).^2 + repmat(FrequencyVector(VolSize(1))',[1 VolSize(2)]).^2),[1 1 VolSize(3)]);
% else
%     RamLak = @(dummy) 1;
% end

% Retrieve refractive index decrement (misnomer 'phase' retrieval on tomogram)
if RamLak
    %hphase = @(a) ifftn( RamLak(size(a)).*PhaseFilter3D( Method, size(a), EnergyDistancePixelsize, RegPar, 0.1, 'single' ).* fftn(a));
    hphase = @(a) ifftn( PhaseFilter3DRamLak( Method, size(a), EnergyDistancePixelsize, RegPar, 0.1, 'single' ).* fftn(a));
else
    hphase = @(a) ifftn( PhaseFilter3D      ( Method, size(a), EnergyDistancePixelsize, RegPar, 0.1, 'single' ).* fftn(a));
end

delta = hcrop( hphase( hpad( vol ) ) );


fprintf(' Phase retrieved using ''%s'' approach in %g min\n',Method,toc/60)
