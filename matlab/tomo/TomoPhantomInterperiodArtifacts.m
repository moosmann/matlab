% Tomographic reconstruction of Shepp-Logan phantom using filtered
% backprojection


%% Phantom and projections
N = 200;
M = 1.5 * N;
NumProj = 1*N;
theta = 0:180/NumProj:180-1/NumProj;
P = phantom('Modified Shepp-Logan',N);
[R, Xp] = radon(P,theta, M);
R = R( round((M-N)/2) + (1:N), :);
RotAxis = ceil(size(R,1)/2);

%% Padding
R = padarray( R, [ceil( 1 * size(R,1)/2), 0], 0, 'both');

%% Inverse radon transformation
freqScal = 1;
OutputSize = size(R,1);
[I1,h] = iradon(R,theta,'linear','Ram-Lak',freqScal,OutputSize);



%% Plotting
subplot(2,2,1), imshow(P), title('Phantom')
subplot(2,2,2), imshow(R, []), title('sino')
subplot(2,2,3), imshow(I1,[]), title('reco')
% 
% %% Phantom and projections
% P = 1-P;
% R = radon(P,theta);
% %% Inverse radon transformation
% freqScal = 1;
% OutputSize = size(R,1);
% I1 = iradon(R,theta,'linear','Ram-Lak',freqScal,OutputSize);
% I2 = iradon(R,theta,'linear','none',freqScal,OutputSize);
% %% Plotting
% figure('Name','Sinogram')
% imshow(R,[])
% figure('Name','Phantom: Original, Filtered BP, Unfiltered BP')
% subplot(1,3,1), imshow(P), title('Original')
% subplot(1,3,2), imshow(I1), title('Filtered backprojection')
% subplot(1,3,3), imshow(I2,[]), title('Unfiltered backprojection')
% 
% %% Writing
% % OutputFolder = '/home/moosmann/data/sim_shepplogan';
% % edfwrite([OutputFolder '/sino/sino_shepplogan_pixel185proj256.edf'],R','float32');
% % imwrite(normat(R),[OutputFolder '/sino/sino_shepplogan_pixel185proj256.tif'])
% % % for ii = 100:-1:1
% % %     R2(:,:,ii) = R;
% % % end
% % for ii=1:size(R,2)
% %     edfwrite(sprintf('%s/proj/proj_%04u.edf',OutputFolder,ii),squeeze(R2(:,ii,:)),'float32');
% % end
% 
% % Check adjoint
% x = ones(444,444);
% NumProj=777; 
% theta = 0:180/NumProj:180-1/NumProj; 
% Ax = radon(x,theta);
% y = ones(size(Ax));
% Ady  = iradon(y,theta,'linear','none',1,max(size(x)));
% Ady = Ady * NumProj/pi * 2;
% 
% n1 = sum(Ax(:) .* y(:));
% n2 = sum(Ady(:) .* x(:));
% fprintf('\n <A x,y>  = %g',n1)
% fprintf('\n <x,Ad y> = %g',n2)
% fprintf('\n <A x,y> / <x,Ad y> = %g',n1/n2-1)
% 
% ishow(Ady)
% ishow(Ax)
% 
