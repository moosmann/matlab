function ScriptTomoBackgroundRemoval()
VolPath = '/export/scratch1/moosmann/ESRF_MI1079_ID19_July2011_inlineTomo/vol/Xenopus_4cell_20keV';

%% Load ART reconstruction with Folkert's background compensation
if 0
    load('/export/scratch1/bleichro/Xenopus_4cell_noring/volume_slice_1001_1100.mat')
    volume = single(volume);
    x = round(size(volume,1)/2) + (-1962/2:1962/2-1);
    volume = volume(x,x,:);
    intf = volume;
    clear volume;
    filename = sprintf('%s/int_folkert_slice1001to1100_%04ux%04ux%04u.mat',VolPath,size(intf));
    save(filename,'intf')
end

%% Load FBP reconstruction without background compensation
if 0
    load(sprintf('%s/fbp_int_zi0925_zf1124_1962x1962x0200.mat',VolPath))
    int0 = vol;clear vol
    int0 = int0(:,:,1001-924 + (0:99));
    filename = sprintf('%s/int_fbp_slice1001to1100_%04ux%04ux%04u.mat',VolPath,size(int0));
    save(filename,'int0')
end

%% Load intensity sino
if 1
   for nn = 100:-1:1
       sino(:,:,nn) = getSino('xeno4cell',12,nn+1000);
   end
   for nn = 100:-1:1
       s = squeeze(sino(:,:,nn));
       s = s - mean(s(:));
       vol(:,:,nn) = iradon(s',360/1600,'linear','Ram-Lak',1,size(sino,2));
   end
   filename = sprintf('%s/int_fbp_sinoMeanSub_slice1001to1100_%04ux%04ux%04u.mat',VolPath,size(vol));
   save(filename,'vol')
   
end



%% Retrieve real refractive index decrement on both tomograms

RetrieveAndSavePhaseMaps(vol,'fbp_sinoMeanSub')
RetrieveAndSavePhaseMaps(intf,'folkert')
RetrieveAndSavePhaseMaps(int0,'fbp')


function RetrieveAndSavePhaseMaps(vol,Name,VolPath)
% padarray along 3rd dimension to avoid jumps of data at boundary which
% introduces ghost arifacts: first and last slices are superimposed
hpad = @(a) padarray(a,[0,0,size(a,3)/2],'symmetric','both');
% crop array back to original dimension
hcrop = @(a) a(:,:,size(a,3)/4+(1:size(a,3)/2));
% loop over different phase retrieval algorithms    
for mm = 1:3
    phaMeth = {'tie','quasi','quasinew'};
    phaMeth = phaMeth{mm};disp(phaMeth)
    % Retrieva refractive index decrement (misnomer 'phase' retrieval on tomogram)
    hphase      = @(a) ifftn( PhaseFilter3D( phaMeth, size(a), [20,0.945,0.75e-6], 2.5, 0.1, 'single' ).* fftn(a));    
    delta = hcrop( hphase( hpad( vol ) ) );
    % save phase map volume
    deltaPath = MakePath('%s/%s_%s_slice1001to1100_%04ux%04ux%04u',VolPath,phaMeth,Name,size(delta));
    for nn = 1:size(delta,3)
        filename = sprintf('%sslice_%04u',deltaPath,nn);
        WriteImage(filename,delta(:,:,nn),'tif')
    end
end