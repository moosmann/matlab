phex = pha_6p00;
tie  = bro_6p00_E10_z0p50_c08000;
ctf  = ctf_6p00_E10_z0p50_c08000;

out=Reco(int_6p00_E10_z0p50_c08000,12,[10 .5 1.1e-6],0.5,1,1);
ctfproj = out.ctfProjected;

x=450:600;

ertie     = abs(phex(x,x)-tie(x,x,1));
ertiepnlo = abs(phex(x,x)-sum(tie(x,x,:),3));
erctf     = abs(phex(x,x)-ctf(x,x,1));
erctfproj = abs(phex(x,x)-ctfproj(x,x));

domain(erctf)
domain(erctfproj)
domain(ertie)
domain(ertiepnlo)

plotRange = [0 2.8];
figure
imshow(ertie,plotRange)
cm = flipud(colormap('gray'));
colorbar;
imtool(erctf,plotRange,'InitialMagnification','fit','Colormap',cm);
imtool(erctfproj,plotRange,'InitialMagnification','fit','Colormap',cm);
imtool(ertie,plotRange,'InitialMagnification','fit','Colormap',cm);
imtool(ertiepnlo,plotRange,'InitialMagnification','fit','Colormap',cm);
