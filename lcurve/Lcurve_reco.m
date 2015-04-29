function Lcurve_reco(data,alphamax,pts,show_plots,filename,padding)
    
    if (nargin<4)||isempty(show_plots)
        show_plots=1;
    end;
    if (nargin<5)||isempty(filename)
        filename='phantom';
    end;
    if (nargin<6)||isempty(padding)
        padding=1;
    end;
    
% Initial alpha.
p0 = 1;
if (p0==0)
    i0=1;
else
    i0=0;
end;
    
% Fourier space cooridnates.
[dimx,dimy] = size(data);
[xi,eta]    = meshgrid(-1/2:1/(dimx):1/2-1/(dimx),-1/2:1/(dimy):1/2-1/(dimy));
xi          = fftshift(xi);
eta         = fftshift(eta);
lap         = xi.^2 + eta.^2;

% Open file.

fid = fopen(sprintf('NumAna_%s.txt',filename), 'w');
% Compute tables of the retrieved phase and correction before and after
% inversion of the Laplacian in dependence of the regularization parameter.
fprintf(fid,['alpha|        <phi1>      var(phi1)|       <phi21>     var(phi21)|       <phi22>' ...
           '     var(phi22)|       <phi23>     var(phi23)|<phi21+phi22>var(phi21+phi22)|          <phi>       var(phi)|\n']);
for ii=pts:-1:p0;
  alpha(ii+i0)         = ii/pts*alphamax;
  [phi_stack,fts,lphi] = rec(data,padding,alpha(ii+i0));
  lphi_stack(:,:,:,ii) = lphi;
  phi1            = phi_stack(:,:,1);
  phi21           = phi_stack(:,:,2);
  phi22           = phi_stack(:,:,3);
  phi23           = phi_stack(:,:,4);
  phi2            = phi21+phi22+phi23;  
  phi             = phi1+phi2;
  fprintf(fid,'%5.4g|%14.8g %14.8g|%14.8g %14.8g|%14.8g %14.8g|%14.8g %14.8g|%14.8g %14.8g|%14.8g %14.8g\n', ...
  alpha(ii+i0), ...
  mean(phi1(:)), max(phi1(:)) -min(phi1(:)), ...
  mean(phi21(:)),max(phi21(:))-min(phi21(:)), ...
  mean(phi22(:)),max(phi22(:))-min(phi22(:)), ...
  mean(phi23(:)),max(phi23(:))-min(phi23(:)), ...
  mean(phi21(:)+phi22(:)),max(phi21(:)+phi22(:))-min(phi21(:)+phi22(:)), ...
  mean(phi(:)),  max(phi(:))  -min(phi(:)));
  
  sol_norm1(ii+i0)   = log10(norm(phi1));
  res_norm1(ii+i0)   = log10(norm(ifft2(lap.*fft2(phi1))-data));
  sol_norm1ft(ii+i0) = log10(norm(fts(:,:,2)));
  res_norm1ft(ii+i0) = log10(norm(lap.*fts(:,:,2)-fts(:,:,1)));
  sol_norm(ii+i0)    = log10(norm(phi));
  res_norm(ii+i0)    = log10(norm(ifft2(lap.*fft2(phi))-data));
  sol_mean1(ii+i0)   = mean(mean(abs(normat(phi1))));
  sol_mean(ii+i0)    = mean(mean(abs(normat(phi))));
  sol_mean1cor(ii+i0)= 10^-alpha(ii)*mean(mean(abs(phi1)));
  sol_meancor(ii+i0) = 10^-alpha(ii)*mean(mean(abs(phi)));

end;
clear ii;

fprintf(fid,['\nalpha|    <L(phi21)>  var(L(phi21))|    <L(phi22)>  var(L(phi22))|' ...
           '<L(phi21+phi22)> var(L(phi21+phi22))\n']);
for ii=p0:pts
fprintf(fid,'%5.3g|%14.8g %14.8g|%14.8g %14.8g|%14.8g %14.8g\n', ...
  alpha(ii+i0), ...
  mean(mean(lphi_stack(:,:,1,ii))),max(max(lphi_stack(:,:,1,ii)))-min(min(lphi_stack(:,:,1,ii))), ...
  mean(mean(lphi_stack(:,:,2,ii))),max(max(lphi_stack(:,:,2,ii)))-min(min(lphi_stack(:,:,2,ii))), ...
  mean(mean(lphi_stack(:,:,1,ii)+lphi_stack(:,:,2,ii))), ... 
  max(max(lphi_stack(:,:,1,ii)+lphi_stack(:,:,2,ii)))-min(min(lphi_stack(:,:,1,ii)+lphi_stack(:,:,2,ii))));    
end;
fprintf(fid,'\n');

if fid~=1,fclose(fid);end;
% L-curve plots.
if show_plots
figure('Name','L-curve: Bronnikov'),plot(res_norm1,sol_norm1,'.-'), ...
for ii=p0:pts, text(res_norm1(ii+i0),sol_norm1(ii+i0),num2str(ii/pts*alphamax));end;
figure('Name','L-curve: Bronnikov + Corrections'),plot(res_norm,sol_norm,'+-'); 
for ii=p0:pts,text(res_norm(ii+i0),sol_norm(ii+i0),num2str(ii/pts*alphamax));end;
figure('Name','L-curve: Bronnikov_ft'),plot(res_norm1ft,sol_norm1ft,'.-');
for ii=p0:pts, text(res_norm1ft(ii+i0),sol_norm1ft(ii+i0),num2str(ii/pts*alphamax));end;
% Norm vs alpha plots.
figure('Name','LogLogPlot of Renormalized Mean vs Alpha: blue=BRO, red=BROCOR'), ... 
  plot(alpha,log10(sol_mean1),'blue',alpha,log10(sol_mean),'red');
figure('Name','Renormalized Mean vs Alpha: blue=BRO, red=BROCOR'), ... 
  plot(alpha,sol_mean1,'blue',alpha,sol_mean,'red');
figure('Name','Regpar corrected Mean vs Alpha: blue=BRO, red=BROCOR'), ... 
  plot(alpha,sol_mean1cor,'blue',alpha,sol_meancor,'red');
end;

