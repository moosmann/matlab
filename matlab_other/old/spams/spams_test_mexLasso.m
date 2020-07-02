clear all
pausetime = 1;
% Data are generated
%X=randn(100,100000);
test_pattern_lena;
lena = SubtractMean( lena );
imsc(lena),
% Add noise
cts=40;
X = 1/cts*double(imnoise(uint16(cts*normat(lena)),'poisson'));
X = SubtractMean(X');
imsc(X),axis tight, pause(pausetime)
% normalise to unit l2-norm
X=X./repmat(sqrt(sum(X.^2)),[size(X,1) 1]);
imsc(X),axis tight, pause(pausetime)
D=lena(:,1:1:end);
imsc(D),axis tight, pause(pausetime)
D=D./repmat(sqrt(sum(D.^2)),[size(D,1) 1]);
imsc(D),axis tight, pause(pausetime)



% parameter of the optimization procedure are chosen
%param.L=20; % not more than 20 non-zeros coefficients (default: min(size(D,1),size(D,2)))
param.lambda=0.15; % not more than 20 non-zeros coefficients
param.numThreads=-1; % number of processors/cores to use; the default choice is -1
                     % and uses all the cores of the machine
param.mode=2;        % penalized formulation

tic
alpha=mexLasso(X,D,param);
t = toc;
fprintf('%f signals processed per second\n',size(X,2)/t);
imsc(D*alpha),axis tight, pause(pausetime)
