InputDir = '/mnt/tomoraid-LSDF/tomo/APS_2BM_LifeCellImaging_GUP28266/lifecellimaging/preprocess/stage11/tomo03_corr';
OutputDir = '/mnt/tomoraid-LSDF/tomo/APS_2BM_LifeCellImaging_GUP28266/lifecellimaging/preprocess/stage11/tomo03_corr_edf';

pha = double(imread(sprintf('%s/phase_0001.tif',InputDir)));
[dimx dimy] = size(pha);
w = FilterGaussian(150,[dimx dimy]);
for nn = 1:1200
    pha = double(imread(sprintf('%s/phase_%04u.tif',InputDir,nn)));
    pha = real(ifft2(w.*fft2(pha)))';
    edfwrite(sprintf('%s/phase_%04u.edf',OutputDir,nn),pha,'float32'); 
end
