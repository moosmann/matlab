function [out,alpha,k,res_norm1,sol_norm1]=Lcurve_2(data,alphamin,alphamax,pts,EnergyDistancePizelsize,show_plots,computePNLO,unprojectedCTF,computeTIElo,BinaryFilter)
% L-curve method to determine regularization parameter.

if nargin<6
    show_plots=1;
end;
if nargin<7
    computePNLO=0;
end
if nargin<8
    unprojectedCTF=0;
end;
if nargin<9
    computeTIElo=1;
end
if nargin<10
    BinaryFilter = 0;
end

padding = 0;

% Fourier space cooridnates.
[dim1,dim2] = size(data);
if padding >0,
    dimx        = padding*2^nextpow2(dim1);
    dimy        = padding*2^nextpow2(dim2);
    data        = padarray(data,[(dimx-dim1)/2,(dimy-dim2)/2],'symmetric','both');
else
    dimx    = dim1;
    dimy    = dim2;
end;
[xi,eta]    = meshgrid(-1/2:1/(dimy):1/2-1/(dimy),-1/2:1/(dimx):1/2-1/(dimx));
xi          = fftshift(xi);
eta         = fftshift(eta);
lap         = xi.^2 + eta.^2;
% Preallocation of memory.
alpha = zeros(1,pts);
sol_norm1 = zeros(1,pts);
res_norm1 = zeros(1,pts);
if computePNLO
    sol_norm2 = zeros(1,pts);
    res_norm2 = zeros(1,pts);
    sol_norm = zeros(1,pts);
    res_norm = zeros(1,pts);
end
tic;
for ii=pts:-1:1;
  alpha(ii) = alphamin + (alphamax - alphamin)*(ii-1)/(pts-1);
  if ii == pts
      out(pts) = RecoGPU(data,alpha(ii),EnergyDistancePizelsize,BinaryFilter,computePNLO,unprojectedCTF,computeTIElo);
  else
      out(ii) = RecoGPU(data,alpha(ii),EnergyDistancePizelsize,BinaryFilter,computePNLO,unprojectedCTF,computeTIElo);
  end
  if computeTIElo
      sol_norm1(ii) = log10(norm(out(ii).tie));
      res_norm1(ii) = log10(norm(ifft2(lap.*fft2(out(ii).tie))-data));
  elseif unprojectedCTF
      sol_norm1(ii) = log10(norm(out(ii).ctf));
      res_norm1(ii) = log10(norm(ifft2(lap.*fft2(out(ii).ctf))-data));          
  elseif BinaryFilter
      sol_norm1(ii) = log10(norm(out(ii).ctfProjected));
      res_norm1(ii) = log10(norm(ifft2(lap.*fft2(out(ii).ctfProjected))-data));      
  end
  if computePNLO
      sol_norm2(ii) = log10(norm(out(ii).tiePNLO));
      res_norm2(ii) = log10(norm(ifft2(lap.*fft2(out(ii).tiePNLO))-data));
      sol_norm(ii)  = log10(norm(out(ii).tie+out(ii).tiePNLO));
      res_norm(ii)  = log10(norm(ifft2(lap.*fft2(out(ii).tie+out(ii).tiePNLO))-data));
  end;
end;
trec = toc;

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
if nargout > 2
    sn2_d1 = diff(abs(sol_norm2),1);
    rn2_d1 = diff(abs(res_norm2),1);
    sn2_d2 = diff(abs(sol_norm2),2);
    rn2_d2 = diff(abs(res_norm2),2);
    sn2_d1 = interp1(1.5:pts,sn2_d1,2:pts-1);
    rn2_d1 = interp1(1.5:pts,rn2_d1,2:pts-1);
    k2     = abs((rn2_d1.*sn2_d2 - rn2_d2.*sn2_d1)./(sn2_d1.^2 + rn2_d1.^2).^1.5);
    i_k2max = find(k2==max(k2));
end
tcur   = toc;
% Print parameters.
fprintf('Time for reconstruction: %gs, time for curvature computation: %gs\n',trec,tcur);
if computePNLO
    fprintf(['At index %i (%i) where alpha = %g (%g) the curvature is maximum ' ...
         'with k = %g (%g)\n'],i_kmax,i_k2max,alpha(i_kmax),alpha(i_k2max),k(i_kmax),k2(i_k2max));
else
       fprintf(['At index %i where alpha = %g the curvature is maximum ' ...
         'with k = %g\n'],i_kmax,alpha(i_kmax),k(i_kmax));
end
% L-curve plots and phase maps.
if show_plots
% Lcurve: Leading order.
figure('Name','L-curve: LO'),plot(res_norm1,sol_norm1,'.-'), ...
for ii=1:pts
    text(res_norm1(ii),sol_norm1(ii),['  (' num2str(ii) ', ' ...
                    num2str(alpha(ii)) ', ' num2str(10^-alpha(ii)) ')']);
end
% Lcurve: Next-to-leading order.
if computePNLO
    figure('Name','L-curve: NLO'),plot(res_norm2,sol_norm2,'.-'), ...
        for ii=1:pts
        text(res_norm2(ii),sol_norm2(ii),['  (' num2str(ii) ', ' ...
            num2str(alpha(ii)) ', ' num2str(10^-alpha(ii)) ')']);
        end
        % Lcurve: LO + NLO
    figure('Name','L-curve: LO + NLO'),plot(res_norm,sol_norm,'.-'), ...
        for ii=1:pts
        text(res_norm(ii),sol_norm(ii),['  (' num2str(ii) ', ' ...
            num2str(alpha(ii)) ', ' num2str(10^-alpha(ii)) ')']);
        end
end
% Curvature of L-curve (LO) vs alpha
figure('Name','Curvature: LO'),plot(2:pts-1,k);
% Curvature of L-curve (NLO) vs alpha
if computePNLO
    figure('Name','Curvature: PNLO'),plot(2:pts-1,k2);
end
end
