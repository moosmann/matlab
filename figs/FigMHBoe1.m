int1=int_0p01_E10_z0p50_c00000;
int2=int_1p00_E10_z0p50_c00000;
int3=int_6p00_E10_z0p50_c00000;

int1 = laf(int1-mean(int1(:)));
int2 = laf(int2-mean(int2(:)));
int3 = laf(int3-mean(int3(:)));

domain(int1)
domain(int2)
domain(int3)

plotRange = [-20 8];%[0 2.8];
FontSize = 32;
cm = (colormap('gray'));
imshow(int1,plotRange,'InitialMagnification','fit','Colormap',cm),colorbar
set(gca,'FontSize',FontSize);

figure,imshow(int1,plotRange,'InitialMagnification','fit','Colormap',cm);
set(gca,'FontSize',FontSize);
figure,imshow(int2,plotRange,'InitialMagnification','fit','Colormap',cm);
set(gca,'FontSize',FontSize);
figure,imshow(int3,plotRange,'InitialMagnification','fit','Colormap',cm);
set(gca,'FontSize',FontSize);


