% Create constant image. Otherwise different noise level would mix when
% mean and variance of the image were computed.
im=ones(2048);
% Apply Poisson noise to image.
NoiseLevel=8000;
imn=double(imnoise(uint16(NoiseLevel*im),'poisson'));
% Print domains.
domain(im,'(image without noise)')
domain(imn,'(image with noise)*1 ')
domain(imn/10,'(image with noise)/10')
domain(10*imn,'(image with noise)*10')