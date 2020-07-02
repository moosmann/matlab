function [k,alpha,res_norm1,sol_norm1,phi1,phi2]=lcurve(data,alphamin,alphamax,pts,lambda,distance,pixelsize,show_plots)

if (nargin<7)||isempty(show_plots),show_plots=1;end

padvalue    = 0;
padding     = 1 ;

% Fourier space cooridnates.
[dim1,dim2] = size(data);
dimx        = padding*2^nextpow2(dim1);
dimy        = padding*2^nextpow2(dim2);
[xi,eta]    = meshgrid(-1/2:1/(dimy):1/2-1/(dimy),-1/2:1/(dimx):1/2-1/(dimx));
xi          = fftshift(xi);
eta         = fftshift(eta);
xieta       = xi.*eta;
lap         = xi.^2 + eta.^2;
% Pad data.
data        = padarray(data,[(dimx-dim1)/2,(dimy-dim2)/2],padvalue,'both');
% Preallocation of data stacks.
phi1        = zeros(dim1,dim2,pts);
phi2        = zeros(dim1,dim2,pts);

tic;
for ii=pts:-1:1
    alpha(ii)       = alphamin + (alphamax - alphamin)*(ii-1)/(pts-1);
    [phi_stack]     = Reco(data,alpha(ii),lambda,distance,pixelsize);
    phi1(:,:,ii)    = phi_stack(:,:,1);
    phi2(:,:,ii)    = phi_stack(:,:,2);
    
    sol_norm1(ii)   = log10(norm(phi_stack(:,:,1)));
    sol_norm(ii)    = log10(norm(phi_stack(:,:,1)+phi_stack(:,:,2)));
    res_norm1(ii)   = log10(norm(ifft2(lap.*fft2(phi1))-data));
    res_norm(ii)    = log10(norm(ifft2(lap.*fft2(phi))-data));
end
trec = toc;

tic;
% Compute curvature.
sn1_d1 = diff(sol_norm1,1);
rn1_d1 = diff(res_norm1,1);
sn1_d2 = diff(sol_norm1,2);
rn1_d2 = diff(res_norm1,2);
sn1_d1 = interp1(1.5:pts,sn1_d1,2:pts-1);
rn1_d1 = interp1(1.5:pts,rn1_d1,2:pts-1);
tcur   = toc;
k      = (rn1_d1.*sn1_d2 - rn1_d2.*sn1_d1)./(sn1_d1.^2 + rn1_d1.^2).^1.5;
i_kmax = 1+find(k==max(k));

% Print parameters.
fprintf('time for reconstruction: %gs, time for curvature computation: %gs\n',trec,tcur);
fprintf('At index %i where alpha = %g the curvature is maximum with k = %g\n',i_kmax,alpha(i_kmax),k(i_kmax-1));

% L-curve plots.
if show_plots
    figure('Name','L-curve: LO'),plot(res_norm1,sol_norm1,'.-'), ...
        for ii=1:pts, text(res_norm1(ii),sol_norm1(ii),num2str(ii/pts*alphamax));end
    figure('Name','Curvature: LO'),plot(2:pts-1,k);
end

