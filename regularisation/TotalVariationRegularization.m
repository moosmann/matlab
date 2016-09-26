%% Source Code for Total Variation Regularization:
clear all
%--------------------------------------------------------------------------
% Total Variation Regularization Method
%--------------------------------------------------------------------------
N = 32;
% make image
I = ones(N);
I = 0.4*I;
a=4;
x=ceil(N/a):ceil(3*N/a);
a=6;
y=ceil(N/a):ceil(3*N/a);
I(x,x)=0.8;
I(y,y)=1;
% blur
[xi,eta]=meshgrid(1/N*((1:N)-floor(N/2)));
I=real(ifft2(fft2(I)./(1+eta.^2+xi.^2)));
ishow(I);
title('Original Image');
xlabel('No:of pixels');
ylabel('No:of pixels');
N = size(I,1); % Image Size
I_res = reshape(I,N^2,1); % Image Vector
% Computing the PSF
b = 10; % Band
s = 7; % Sigma value
Z = [exp(-((0:b-1).^2)/(2*s^2)),zeros(1,N-b)];
A = toeplitz(Z);
A = (1/(2*pi*s^2))*kron(A,A); % Point Spread Function
% Addition of Gaussina Noise
noiselevel = 0.0001;
noise = imnoise(I,'gaussian',0,noiselevel);
ishow(noise);
title('Noise added to the image');
xlabel('No:of pixels');
ylabel('No:of pixels');
% % Creating Noise Image.
% x = zeros(N,N);
% noisex = imnoise(x,'gaussian',noiselevel);
% ishow(noisex);
% title('Noise');
% Lexicographic Arrangement of Noisy image
Nr = reshape(noise,N^2,1);
lg_Nr = double(Nr);
% g is the blurred image
% I_res is the resized image to N*1
% I_org is the image restored back to original size of N*N
g = A*lg_Nr;
I_blur = reshape(g,N,N);
ishow(I_blur);
title('Blurred and Noisy Image');
xlabel('No:of pixels');
ylabel('No:of pixels');
% Computing the Least square method
% using in-built command
lsinbuilt = lsqr(A,g,1e-06,100);
ishow(reshape(lsinbuilt,N,N));
title('Least Square Output using the inbuilt Matlab command');
xlabel('No:of pixels');
ylabel('No:of pixels');
% Least suqares using formulation
%inverse = inv(transpose(A)*A)*transpose(A)*g;lsrest = reshape(inverse,N,N);
lsrest = reshape(inv(transpose(A)*A)\(transpose(A)*g),N,N);
ishow(lsrest);
title('Least Square Output using the formulation');
xlabel('No:of pixels');
ylabel('No:of pixels');
fixed_iter = 24; % fixed point iteration
beta = 0.1; % Smoothing factor
lamda = 1.5e-05; % Regularization parameter
% Non-Linear Equation
% Computation of Regularization operator:
n = N;
nsq = n^2;
Delta_x = 1 / n;
Delta_y = Delta_x;
D = spdiags([-ones(n,1) ones(n,1)], [0 1], n,n) / Delta_x;
I_trunc1 = spdiags(ones(n,1), 0, n,n);
Dx1 = kron(D,I_trunc1); % Forward (upwind) differencing in x
Dy1 = kron(I_trunc1,D); % Forward (upwind) differencing in y
f_fp = zeros(n,n); % Initial Guess of Blank Image
fvec = f_fp(:); % Image Vector
i=1;
cng(fixed_iter) = 0;
time(fixed_iter) = 0;
for f_iter = 1:fixed_iter
tic
t=((Dx1*fvec).^2 + (Dy1*fvec).^2);
psi_prime1 = 1./sqrt(t + (0.1)^2);
Dpsi_prime1 = spdiags(psi_prime1, 0, (n)^2,(n)^2);
L1 = Dx1' * Dpsi_prime1 * Dx1 + Dy1' * Dpsi_prime1 * Dy1;
L = L1 * Delta_x * Delta_y;
pen_fun = lamda * L * fvec; % Computing alpha*operator*initial image
% computing the gradient of the non-Linear equation
Gv = A' * (A * fvec - g) + pen_fun;
Gv_nm = norm(Gv); % Compute the norm of the fixed point gradient
H_residual = lamda * L; % Computation alpha * operator
H = (A' * A) + H_residual; % Computes the Hessian Matrix
S = -(inv(H)) * Gv; % Compute the quasi-newton
fvec = fvec + S; % Upgrading the initial image
df = double(I_res) - double(fvec);
cng(f_iter) = norm(df);
toc
time(f_iter) = toc;
end
semilogy(1:fixed_iter, Gv_nm,'r-');
title('Gradient plot of the fixed point iteration-Convergence');
xlabel('No:of iteration');
ylabel('||gradient||');
grid on
semilogy(1:fixed_iter,cng,'r-');
title('Difference between the True image and the restored image (Lamda=1.5e-05)');
xlabel('No:of iteration');
ylabel('|| True Image - Restored Image ||');
grid on
y = [2.5882, 2.5837, 2.5825, 2.5816, 2.5808, 2.5801, 2.5795, 2.5789, 2.5785,2.5780, 2.5778, 2.5775, 2.5772, 2.5771, 2.5769, 2.5768, 2.5767,2.5766, 2.5766, 2.5762, 2.5773, 2.5794, 2.5824, 2.5857, 2.5923, 2.5991,2.6054, 2.6116, 2.6180, 2.6242];
x = [10^-6, 1.5*10^-6, 2*10^-6, 2.5*10^-6, 3*10^-6, 3.5*10^-6, 4*10^-6, 4.5*10^-6, 5*10^-6, 5.5*10^-6, 6*10^-6, 6.5*10^-6, 7*10^-6, 7.5*10^-6, 8*10^-6, 8.5*10^-6, 9*10^-6, 9.5*10^-6, 1*10^-5, 1.5*10^-5, 2.0*10^-5, 2.5*10^-5, 3.0*10^-5, 3.5*10^-5, 4.5*10^-5, 5.5*10^-5, 6.5*10^-5, 7.5*10^-5,8.5*10^-5, 9.5*10^-5];
figure, plot(x,y)
xlim([1e-06 4e-05])
xlabel('Regularization Parameter (Lamda)');
ylabel('|| True Image - Restored Image ||');
title('Selection of Optimal Lamda');
grid on
fvec_res = reshape (fvec, N, N);
ishow(fvec_res);
title ('Restored Image: Total Variation process');
xlabel ('No: of pixels');
ylabel ('No: of pixels');