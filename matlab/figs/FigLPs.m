x = 1:1000;
y = 532;
phex = pha_0p15;
ctf1 = ctf_0p15_E10_z0p10_c00000(:,:,1);
ctf2 = ctf_1p00_E10_z0p10_c00000(:,:,1);
ctf3 = ctf_6p00_E10_z0p10_c00000(:,:,1);
%figure,plot(x,phex(x,y),'black',x,ctf1(x,y),'green',x,ctf2(x,y),'blue',x,ctf3(x,y),'red');
figure,plot(x,phex(x,y),'black',x,ctf1(x,y),'green',x,ctf2(x,y),'blue');
