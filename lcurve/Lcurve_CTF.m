function [lo,nlo,k,alpha,res_norm1,sol_norm1]=Lcurve_CTF(data,alphamin,alphamax,pts,cor,lambda,distance,pixelsize,show_plots,write_res)
% L-curve method to determine regularization parameter.

if nargin<5
    cor =1;
end;
if nargin<6
    lambda = 1;
end;
if nargin<7
    distance = 1;
end;
if nargin<8
    pixelsize = 1;
end;
if nargin<9
    show_plots=1;
end;
if nargin<10
    write_res=0;
end;

padding = 0 ;

% Fourier space cooridnates.
[dim1,dim2] = size(data);
if padding >0,
    dimx        = padding*2^nextpow2(dim1);
    dimy        = padding*2^nextpow2(dim2);
    data        = padarray(data,[(dimx-dim1)/2,(dimy-dim2)/2],padvalue,'both');
else
    dimx    = dim1;
    dimy    = dim2;
end;
[xi,eta]    = meshgrid(-1/2:1/(dimy):1/2-1/(dimy),-1/2:1/(dimx):1/2-1/(dimx));
xi          = fftshift(xi);
eta         = fftshift(eta);
%xieta       = xi.*eta;
lap         = xi.^2 + eta.^2;
% Preallocation of data stacks.
lo        = zeros(dimx,dimy,pts);
if cor==1,
    nlo    = zeros(dimx,dimy,pts);
end;

tic;
for ii=pts:-1:1;
  alpha(ii)     = alphamin + (alphamax - alphamin)*(ii-1)/(pts-1);
  [phi_stack]   = RecoCTF(data,alpha(ii),lambda,distance,pixelsize,1,0,0,cor);
  lo(:,:,ii)  = phi_stack(:,:,1);
  if cor==1,
      nlo(:,:,ii)  = phi_stack(:,:,2);
  end;
  
  sol_norm1(ii) = log10(norm(phi_stack(:,:,1)));
  res_norm1(ii) = log10(norm(ifft2(lap.*fft2(phi_stack(:,:,1)))-data));

  sol_norm2(ii) = log10(norm(phi_stack(:,:,2)));
  res_norm2(ii) = log10(norm(ifft2(lap.*fft2(phi_stack(:,:,2)))-data));

  if cor==1,
      sol_norm(ii)  = log10(norm(phi_stack(:,:,1)+phi_stack(:,:,2)));
      res_norm(ii)  = log10(norm(ifft2(lap.*fft2(phi_stack(:,:,1)+phi_stack(:,:,2)))-data));
  end;
end;
trec = toc;
if cor==0,
    nlo = lo;
end;

tic;
% Compute curvature.
sn1_d1 = diff(abs(sol_norm1),1);
rn1_d1 = diff(abs(res_norm1),1);
sn1_d2 = diff(abs(sol_norm1),2);
rn1_d2 = diff(abs(res_norm1),2);
sn1_d1 = interp1(1.5:pts,sn1_d1,2:pts-1);
rn1_d1 = interp1(1.5:pts,rn1_d1,2:pts-1);
k      = abs((rn1_d1.*sn1_d2 - rn1_d2.*sn1_d1)./(sn1_d1.^2 + rn1_d1.^2).^1.5);
i_kmax = find(k==max(k));
if cor==1,
    sn2_d1 = diff(abs(sol_norm2),1);
    rn2_d1 = diff(abs(res_norm2),1);
    sn2_d2 = diff(abs(sol_norm2),2);
    rn2_d2 = diff(abs(res_norm2),2);
    sn2_d1 = interp1(1.5:pts,sn2_d1,2:pts-1);
    rn2_d1 = interp1(1.5:pts,rn2_d1,2:pts-1);
    k2     = abs((rn2_d1.*sn2_d2 - rn2_d2.*sn2_d1)./(sn2_d1.^2 + rn2_d1.^2).^1.5);
    i_k2max = find(k2==max(k2));
end;
tcur   = toc;

% Print parameters.
fprintf('Time for reconstruction: %gs, time for curvature computation: %gs\n',trec,tcur);
if cor==1,
    fprintf(['At index %i (%i) where alpha = %g (%g) the curvature is maximum ' ...
         'with k = %g (%g)\n'],i_kmax,i_k2max,alpha(i_kmax),alpha(i_k2max),k(i_kmax),k2(i_k2max));
else
       fprintf(['At index %i where alpha = %g the curvature is maximum ' ...
         'with k = %g\n'],i_kmax,alpha(i_kmax),k(i_kmax));
end;

% L-curve plots and phase maps.
if show_plots
% Lcurve: Leading order.
figure('Name','L-curve: LO'),plot(res_norm1,sol_norm1,'.-'), ...
for ii=1:pts, text(res_norm1(ii),sol_norm1(ii),['  (' num2str(ii) ', ' ...
                    num2str(alpha(ii)) ', ' num2str(10^-alpha(ii)) ')']);end;
if write_res,saveas(gcf,[folder 'Lcurve_LO.tif'],'tiffn');end;
% Lcurve: Next-to-leading order.
if cor==1,
figure('Name','L-curve: NLO'),plot(res_norm2,sol_norm2,'.-'), ...
for ii=1:pts, text(res_norm2(ii),sol_norm2(ii),['  (' num2str(ii) ', ' ...
                    num2str(alpha(ii)) ', ' num2str(10^-alpha(ii)) ')']);end;
if write_res,saveas(gcf,[folder 'Lcurve_LO.tif'],'tiffn');end;
% Lcurve: LO + NLO
figure('Name','L-curve: LO + NLO'),plot(res_norm,sol_norm,'.-'), ...
for ii=1:pts, text(res_norm(ii),sol_norm(ii),['  (' num2str(ii) ', ' ...
                    num2str(alpha(ii)) ', ' num2str(10^-alpha(ii)) ')']);end;
if write_res,saveas(gcf,[folder 'Lcurve_NLO.tif'],'tiffn');end;
end;
% Curvature of L-curve (LO) vs alpha
figure('Name','Curvature: LO'),plot(2:pts-1,k);
if write_res,saveas(gcf,[folder 'CurvatureOfLcurve_LO.tif'],'tiffn');end;
% Curvature of L-curve (NLO) vs alpha
if cor==1,figure('Name','Curvature: NLO'),plot(2:pts-1,k2);end;

% Phase
figure('Name','Phase map at maximimum curvature of Lcurve'),
imshow(lo(:,:,i_kmax),[],'InitialMagnification','fit');
if write_res,edfwrite([folder 'PhaseMapAtMaxCurvatureOfLcurve.edf'],lo(:,:,i_kmax),'float32');end
end;

