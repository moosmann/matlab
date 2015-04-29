function [mean_bron,mean_improved] = BinsengrasRec(padding,alpha,beta,fig);
% Phase retrieval for BINSENGRAS. 

global proj_bins;

if isempty(proj_bins)
[head,proj_bins] = pmedf_read('binsengras_projection.edf');
proj_bins        = transpose(proj_bins);
end;

%compute Bronnikov and corrections
[phi0,phi1,phi2] = rec(proj_bins,padding,alpha,beta);
phi    = (phi0 + phi1 + phi2);

%compute mean
mean_bron     = mean(mean(abs(normat(phi0))));
mean_improved = mean(mean(abs(normat( phi))));

%figures
if fig==1
%figure,imshow(proj_bins,[],'InitialMagnification','fit'),colorbar;
figure('Name','Binsengras: phi0(Bronnikov reconstruction)'), ... 
    imshow(phi0,[],'InitialMagnification','fit'),colorbar;
end;

if fig==2

figure('Name','Binsengras: correction phi1'),  ...
    imshow(phi1,[],'InitialMagnification','fit'),colorbar;
figure('Name','Binsengras: correction phi2'), ... 
    imshow(phi2,[],'InitialMagnification','fit'),colorbar;
%figure,imshow(phi1,[],'InitialMagnification','fit'),colorbar;
%figure,imshow(phi,[],'InitialMagnification','fit'),colorbar;
end;

if fig==3
figure('Name','Binsengras: phi0 + phi1 +phi2'), ... 
    imshow(phi,[],'InitialMagnification','fit'),colorbar;
end;

if fig==4
figure,imshow([normat(proj_bins),normat(phi0),normat(phi1)],[], ...
	      'InitialMagnification','fit');
end;

