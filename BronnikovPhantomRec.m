function [mean_bron,mean_improved,diff0,diff] = BronnikovPhantomRec(proj,slice, ...
                                                  padding,alpha,beta,fig);
% Phase retrieval for the Bronnikov phantom to first and second
% order. Computation of mean maps and error maps. Optional printing of
% figures for several (phase) maps.

global phase;

pha           = phase(:,:,slice);

%compute Bronnikov and corrections
[phi0,phi1,phi2] = rec(proj(:,:,slice),padding,alpha,beta);
phi  = phi0 + phi1 + phi2;

%mean and error estimate for bronnikov phantom
mean_bron     = mean(mean(abs(normat(phi0))));
mean_improved = mean(mean(abs(normat(phi ))));
diff0         = mean(mean(abs(normat(pha)-normat(phi0))));
diff          = mean(mean(abs(normat(pha)-normat(phi ))));

if 0
minval = min(min(phi));
maxval = max(max(phi));
phi0   = (phi0-minval)/(maxval-minval);
phi1   = (phi1)/(maxval-minval);
phi2   = (phi2)/(maxval-minval);
ranges(phi);
ranges(phi0);
ranges(phi1);
ranges(phi0+phi1);
ranges(phi2);
display(minval);
display(maxval);
end;

%figures
if fig>0

if fig==1
figure('Name','Phantom: phi0'), ... 
    imshow(phi0,[],'InitialMagnification','fit'),colorbar;
end;

if fig==2
figure('Name','Phantom: first and second correction'),
 imshow(phi1,[],'InitialMagnification','fit'),colorbar;
figure('Name','Phantom: third correction'),
 imshow(phi2,[],'InitialMagnification','fit'),colorbar;
figure('Name','Phantom: second + third correction'),
 imshow(phi1+phi2,[],'InitialMagnification','fit'),colorbar;
end;

if fig==3
figure('Name','Phantom: phase: bronnikov'), ...
       imshow([phi0],[],'InitialMagnification','fit'),colorbar;  
figure('Name','Phantom: phase: bronnikov + corrections'), ...
       imshow([phi],[],'InitialMagnification','fit'),colorbar;  
end;

if fig==4
phaphi0   = normat(pha)-normat(phi0);
phaphi    = normat(pha)-normat(phi);
figure('Name','Phantom; phase_exact-phase_bron, phase_exact-phase_improved'), ...
       imshow([abs(phaphi0),abs(phaphi)],[],'InitialMagnification','fit'),colorbar;
end;

if fig==5
figure('Name','Phantom: phase: bronnikov and bronnikov plus corrections'), ...
       imshow([phi0,phi],[],'InitialMagnification','fit'),colorbar;  
end;


end;



