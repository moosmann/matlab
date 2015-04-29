im = padarray(0.5+0.5*normat(double(imread('/home/moosmann/lena.tif'))),[256 ...
                   256],'replicate','both');

sigma = 16;
win = padarray(ones(512+8*sigma),[256-4*sigma 256-4*sigma],0,'both');
win = imfilter(win,fspecial('gaussian',[256 256],sigma));

im = win.*im;

edp = [10 .5 1e-6];
padding = 2;
MethodOrPadValue = 0;
PrintDomains = 0;
int = Propagation(im,edp,padding,MethodOrPadValue,PrintDomains);
im = im - mean(im(:));



clear sigma padding MethodOrPadValue PrintDomains

[tie,ctf]=Reco(int,12,edp,0.5,1);

er0=abs(tie(:,:,1)-im);er1=abs(tie(:,:,1)+tie(:,:,2)-im);erpctf=abs(ctf(:,:,2)-im);
