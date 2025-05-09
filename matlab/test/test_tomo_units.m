
%% Phantom and projections
N = 1024;
NumProj = 4*N;
theta = 0:180/NumProj:180-1/NumProj;
P = phantom('Modified Shepp-Logan',N);
[R, Xp] = radon(P,theta);
RotAxis = ceil(size(R,1)/2);
R =  single(R);
write32bitTIFfromSingle("/home/moosmanj/beamtimeid/processed/scan/sino/rawBin1/sino.tif",R');

%% Inverse radon transformation
freqScal = 1;
OutputSize = size(R,1);
[I1,h] = iradon(single(R),theta,'linear','Ram-Lak',freqScal,OutputSize);


%% Plotting
figure('Name','Phantom: Original, FBP, BP')
subplot(1,4,1)
imshow(R,[])
title('sinogram')
colorbar

subplot(1,4,2)
imshow(P,[])
title('Original')
colorbar

subplot(1,4,3)
imshow(I1,[])
title('FBP')
colorbar


subplot(1,4,4)
imshow(I1,[])
title('FBP')
colorbar

