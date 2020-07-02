function [phialpha1,phialpha2] = MeanVsAlpha(data,alphamin,alphamax,pts,padding);
% Plot the mean of the absolute value of the (renormalized) phase map for LO
% and LO+NLO versus the regularisaton parameter alpha. The retrieved phase
% has to be regularized since the absolute value of retrieved phase strongly
% depends on the regularization parameter.

lambda = 1;
distance = 1;
pixelsize = 1;

if (nargin<4),padding=1;end;
padding =1;
% Reconstruction of phase maps of different alphas.
for ii=0:pts;
  x(ii+1)       = alphamin+ii/pts*(alphamax-alphamin);
  phi           = Reco(data,x(ii+1),lambda,distance,pixelsize,padding);
  phi1          = phi(:,:,1);
  phi2          = phi(:,:,2);
  phi           = phi1 + phi2;
  phialpha1(:,:,ii+1) = phi1;
  phialpha2(:,:,ii+1) = phi2;
  meanb(ii+1)   = mean(abs(normat(phi1(:))));
  meanbc(ii+1)  = mean(abs(normat(phi(:))));
  meanb2(ii+1)  = mean(abs(normat(phi1(:))).^2);
  meanbc2(ii+1) = mean(abs(normat(phi(:))).^2);
  meanub(ii+1)   = mean(abs(phi1(:)));
  meanubc(ii+1)  = mean(abs(phi(:)));
  meanub2(ii+1)  = mean(abs(phi1(:)).^2);
  meanubc2(ii+1) = mean(abs(phi(:)).^2);
end;
% Figure: mean0,mean1
int=1:pts+1;
figure('Name','Mean vs Alpha: normed lo, normed lo&nlo; lo, lo&nlo (blue=BRO, red=BROCOR)'),
subplot(2,2,1),plot(x(int),meanb(int),'blue',x(int),meanbc(int),'red'),
subplot(2,2,2),plot(x(int),meanb2(int),'blue',x(int),meanbc2(int),'red'),
subplot(2,2,3),plot(x(int),meanub(int),'blue',x(int),meanubc(int),'red'),
subplot(2,2,4),plot(x(int),meanub2(int),'blue',x(int),meanubc2(int),'red'),

