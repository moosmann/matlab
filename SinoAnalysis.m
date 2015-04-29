
N = 1024;
[x,y] = meshgrid((1:N)-N/2,(1:N)-N/2);
m0 = double((y == round(N/4*sin(2*pi*x/N))));
m1 = double((y == round(N/2*sin(2*pi*x/N))));
m2 = double((y == round(1*N*sin(2*pi*x/N))));
m3 = double((y == round(2*N*sin(2*pi*x/N))));

M = 100;
for nn = M:-1:1
    m(:,:,nn) = double((y == round(10*N/2*(nn/M)*sin(2*pi*x/N))));    
    
    fdir = 1;
    %mf1(:,:,nn) = log(1+abs(fftshift(fft(squeeze(m(:,:,nn)),[],fdir),fdir)));
    %mf1(:,:,nn) = log(1+abs(fft(squeeze(m(:,:,nn)),[],fdir)));
    %mf1(:,:,nn) = (angle(fft(squeeze(m(:,:,nn)),[],fdir)));
    md(:,:,nn) = dst(squeeze(m(:,:,nn))');
    
    fdir = 2;
    %mf2(:,:,nn) = log(1+abs(fftshift(fft(squeeze(m(:,:,nn)),[],fdir),fdir)));
    %mf2(:,:,nn) = (angle(fft(squeeze(m(:,:,nn)),[],fdir)));
    
    %mfb(:,:,nn) = (angle(fftshift(fft2(squeeze(m(:,:,nn))))));
end



%itool(m)
%ftool(m0,0,0,0)