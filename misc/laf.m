function out = laf(int)
    
int = int - mean(int(:));
out = log(abs(fftshift(fft2(int))));
