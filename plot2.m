function plot2

    padding = 1;
    alpha   = 12;
    iterations = 0;
global phase proj40 proj50;

pha_30           = normat(phase(:,:,60));
pha_135          = normat(phase(:,:,270));

[phi0,phi1,phi2] = rec(proj40(:,:,60),padding,alpha,iterations);
phi0_30 = normat(phi0);
phi_30 = normat(phi0+phi1+phi2);
[phi0,phi1,phi2] = rec(proj40(:,:,270),padding,alpha,iterations);
phi0_135 = normat(phi0);
phi_135 = normat(phi0+phi1+phi2);

[phi0,phi1,phi2] = rec(proj50(:,:,60),padding,alpha,iterations);
phi50_0_30 = normat(phi0);
phi50_30 = normat(phi0+phi1+phi2);
[phi0,phi1,phi2] = rec(proj50(:,:,270),padding,alpha,iterations);
phi50_0_135 = normat(phi0);
phi50_135 = normat(phi0+phi1+phi2);

a=max(max(pha_30-phi0_30));
balk=a*ones(256,20);

figure,imshow([...
a*(1-pha_30),balk,pha_30-phi0_30,pha_30-phi_30,balk,pha_30-phi50_0_30,pha_30-phi50_30;...
a*(1-pha_135),balk,pha_135-phi0_135,pha_135-phi_135,balk,pha_135-phi50_0_135,pha_135-phi50_135...
],[],'InitialMagnification','fit')
saveas(gcf,sprintf('phantom_errorplots_d40_d50_th30_th135.eps'),'psc2');