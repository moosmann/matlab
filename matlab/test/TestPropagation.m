%% Test forward propagaton using Fresnel propagator
fprintf('\nStart\n')
ca
edp = [10 .3 2e-6];
sf = 1/1000;
%% Load lena test pattern
dim=1024;
dimX = floor(dim/2);
dimY = dimX;
[xi eta] = meshgrid(0:dimX-1,0:dimY-1);
testmap = abs(1/2+cos(pi*xi/(dimX-1))/2).*abs(1/2+cos(pi*eta/(dimY-1))/2);
testmap = padarray(testmap,[round((dim-dimX)/2) round((dim-dimY)/2)],0,'both');
%testmap = phantom('Shepp-Logan');%LenaMap;
%% Propagation
[intPurePhase] = Propagation2(testmap,0,edp,1);
[intPureAbs]   = Propagation2(0,testmap,edp,1);
[intMixed]   = Propagation2( testmap,sf*testmap,edp,1);
%% Show
domain(testmap);
ishow(testmap)
domain(intPurePhase);
ishow(intPurePhase);
%domain(intPureAbs);
%ishow(intPureAbs);
domain(intMixed);
ishow(intMixed);
% x=dimX/2+(-20:20);
% y=dimX;
% figure,
% plot(x,intPurePhase(x,y),x,intPureAbs(x,y));
% fprintf('\nFinished')