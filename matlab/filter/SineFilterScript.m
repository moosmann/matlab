function [phi,int,f,intFTfilt] = SineFilterScript(int,Ezp,height)

   if nargin<3,height=0.1;end;

show = 0;

resolution = size(int);
f          = SineFilter(resolution,Ezp,height);
if show, ishow(f);ishow(int); end;
intFTfilt  = fftshift(fft2(int)).*f;
int        = ifft2(fftshift(intFTfilt));
if show, ishow(int); end;
phi        = RecoCTF(int,12,Ezp,0);
if show, ishow(phi,2); end;
