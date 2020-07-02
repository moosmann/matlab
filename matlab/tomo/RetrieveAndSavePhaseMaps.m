function delta = RetrieveAndSavePhaseMaps(vol,Name,VolPath,RegPar,doWrite)

if nargin < 4
    RegPar = 2.5;
end
if nargin < 5
    doWrite = 0;
end

% padarray along 3rd dimension to avoid jumps of data at boundary which
% introduces ghost arifacts: first and last slices are superimposed
hpad = @(a) padarray(a,[0,0,size(a,3)/2],'symmetric','both');
% crop array back to original dimension
hcrop = @(a) a(:,:,size(a,3)/4+(1:size(a,3)/2));
% loop over different phase retrieval algorithms    
for mm = 2
    phaMeth = {'tie','quasi','quasinew'};
    phaMeth = phaMeth{mm};disp(phaMeth)
    % Retrieva refractive index decrement (misnomer 'phase' retrieval on tomogram)
    hphase      = @(a) ifftn( PhaseFilter3D( phaMeth, size(a), [20,0.945,0.75e-6], RegPar, 0.1, 'single' ).* fftn(a));
    delta = hcrop( hphase( hpad( vol ) ) );
    % save phase map volume
    if doWrite > 0
        deltaPath = MakePath('%s/%s_%s_slice1001to1100_%04ux%04ux%04u',VolPath,phaMeth,Name,size(delta));
        for nn = 1:size(delta,3)
            filename = sprintf('%sslice_%04u',deltaPath,nn);
            WriteImage(filename,delta(:,:,nn),'tif')
        end
    end
end