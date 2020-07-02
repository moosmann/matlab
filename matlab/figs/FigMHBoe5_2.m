phex = pha_6p00;
tie  = bro_6p00_E10_z0p50_c08000;
ctf  = ctf_6p00_E10_z0p50_c08000;

x=450:600;

ertie     = abs(phex(x,x)-tie(x,x,1));
ertiepnlo = abs(phex(x,x)-sum(tie(x,x,:),3));
erctf     = abs(phex(x,x)-ctf(x,x,1));
erctfproj = abs(phex(x,x)-ctfproj(x,x));

domain(ertie)
domain(ertiepnlo)
domain(erctf)
domain(erctfproj)

cm = (colormap('gray'));
imtool(ertie,[],'InitialMagnification','fit','Colormap',cm);
imtool(ertiepnlo,[],'InitialMagnification','fit','Colormap',cm);
