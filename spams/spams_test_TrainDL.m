%clear all; 
patchSize = [8 8];

% Read images
I = double( imread( [ getenv('HOME') '/data/test_pattern/lena/lena.png' ] ) )/255;
imSize = size( I );
I2 = double( imread( [ getenv('HOME') '/data/test_pattern/barbara/barbara_512x512.png' ] ) )/255;
I3 = double( imread( [ getenv('HOME') '/data/gate/vol/20150320_shifted_water_spheres/vol_z0050.tif' ] ) );
I20 = I2;
cts = 100;
I2 = 1 / cts * double( imnoise( uint16( cts * ( I20 ) ), 'poisson' ) );

% extract 8 x 8 patches
X = im2col(I,   patchSize,'sliding');
X2 = im2col(I2, patchSize,'sliding');
X3 = im2col(I3, patchSize,'sliding');

%X=X-repmat(mean(X),[size(X,1) 1]);
%X=X ./ repmat(sqrt(sum(X.^2)),[size(X,1) 1]);

param.K=100;  % learns a dictionary with K elements
param.lambda=0.8;
param.numThreads=-1; % number of threads
param.batchsize=400;
param.verbose=false;
param.iter=1000;  % let us see what happens after 1000 iterations.

%% Train dictionary
fprintf('Dictionary Learning. ');
tic
D = mexTrainDL(X,param);
%D = mexTrainDL_Memory(X,param);
fprintf('Elapsed time: %f\n',toc);

% Show
ImD = displayPatches(D);
subplot(4,2,1);
imagesc(ImD); colormap('gray'); axis equal tight
title('Dictionary patches')

%% Find sparse representation
%param.lambda=0.15;
fprintf('Evaluating cost function: lena. \n');
alpha = mexLasso( X, D, param);
R = mean(0.5*sum((X-D*alpha).^2)+param.lambda*sum(abs(alpha)));
fprintf('Elapsed time: %f\n',toc);

fprintf('Evaluating cost function: barbara. \n');
alpha2 = mexLasso( X2, D, param);
R2 = mean(0.5*sum((X2-D*alpha2).^2)+param.lambda*sum(abs(alpha2)));
fprintf('Elapsed time: %f\n',toc);

fprintf('Evaluating cost function: gate spheres. \n');
alpha3 = mexLasso( X3, D, param);
R3 = mean(0.5*sum((X3-D*alpha3).^2)+param.lambda*sum(abs(alpha3)));
fprintf('Elapsed time: %f\n',toc);

fprintf('Objective functions: %f (lena), %f (barbara)\n', R, R2);

%% Sparse representation alpha in dictionary D
fprintf('Rearrange images. \n');
xalpha = D * alpha;
xalpha2 = D * alpha2;
xalpha3 = D * alpha3;

% Rearrange matrix into image
Imxalpha = spams_col2im( xalpha, patchSize, imSize );
Imxalpha2 = spams_col2im( xalpha2, patchSize, imSize );
Imxalpha3 = spams_col2im( xalpha3, patchSize, size(I3) );
fprintf('Elapsed time: %f\n',toc);

% Show
subplot(4,2,2);
imagesc(Imxalpha); colormap('gray'); axis equal tight
title('Sparse repr of training data lena')
subplot(4,2,3);
imagesc(I2); colormap('gray'); axis equal tight off
title('Noisy data barbara')
subplot(4,2,4);
imagesc(Imxalpha2); colormap('gray'); axis equal tight off
title('Sparse repr of noisy data')
subplot(4,2,5);
imagesc(I2-I20); colormap('gray'); axis equal tight off
title('Noisy data - ground truth')
subplot(4,2,6);
imagesc(Imxalpha2-I20); colormap('gray'); axis equal tight off
title('Sparse repr of noisy data - ground truth')

subplot(4,2,7);
imagesc(I3); colormap('gray'); axis equal tight off
title('gate shpere')
subplot(4,2,8);
imagesc(Imxalpha3); colormap('gray'); axis equal tight off
title('Sparse repr of gate sphere')

drawnow;