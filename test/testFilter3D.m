% Test phase reconstruction on 3D volume. Phase and intensity maps were
% simulated by Daniel using SRCL. The projection were reconstructed using
% PyHST.

% Prefix: ice_cone_
% range: 0-179.5 degree in 0.5 degree steps
% propagation distance: 0.4m
% energy: 30 keV
% pixel size: 1.4um
% material: ice
% rotation axis position x=127 (besser prüfen, aber die Indizierung sollte von 0 bis 255 laufen)
% resolution: 256 x 256


ProjPath = '/mnt/ANKA-LSDF-DDN/tomoraid3/user/haenschke/simulations/volPhaseRetrieval/';
if ~exist('vol','var');
    vol = mexVolRead([ProjPath 'vol/int.vol'],[256 256 256],'float32');
end


out = Filter3D(vol,2.5,[30 0.4 1.4e-6]);

