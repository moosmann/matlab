function vol = iradonLoop(sino,VolPath,Name)


for nn = size(sino,3):-1:1
    s = squeeze(sino(:,:,nn));
    s = padarray(s,[0 size(s,2)/2],'replicate','both');
    s = s - mean(s(:));
    vol(:,:,nn) = iradon(s',360/1600,'linear','Ram-Lak',1,size(sino,2));
    % Write images
    deltaPath = MakePath('%s/int_%s_slice1001to1100_%04ux%04ux%04u',VolPath,Name,size(vol));    
    filename = sprintf('%sslice_%04u',deltaPath,nn);
    WriteImage(filename,vol(:,:,nn),'tif')
    
end
