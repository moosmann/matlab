function [intstack,phi1stack,phi2stack,err1stack,err2stack,merr]= ...
        Distance(phase_object,z_max,z_pts,energy,pixelsize,alpha,padding,normpow,errfun)
    
    if(nargin<7),padding=1;end;
    if(nargin<8),normpow=1;end;
    if(nargin<9),errfun=1;end;

padvalue    = 0;
renormalize = 0;
lambda      = EnergyConverter(energy);
% Print parameters.
fprintf(1,'distance_max=%g, lambda=%g, pixelsize=%g, lambda*distance_max/pixelsize^2=%g\n', ...
        z_max,lambda,pixelsize,lambda*z_max/pixelsize^2);
% Print domains.
phase_object = phase_object - mean(phase_object(:));
pmax         = max(phase_object(:));
pmin         = min(phase_object(:));
fprintf(1,['Domain of phase object:[%1.9g,%1.9g], Mean=%1.9g, Max-Min=%1.9g\n'], ...
        pmin,pmax,mean(phase_object(:)),pmax-pmin);
% Loop over distances.
p0=1;
if (p0==0),i0=1;else,i0=0;end;
for ii=p0:z_pts,
    dist     = ii/z_pts*z_max;
    z(ii+i0) = dist;
    [int]    = Propagation(phase_object,dist,lambda,pixelsize,padding,0);
    [phi]    = Reco(int,alpha,padding,padvalue,renormalize);
    phi1     = phi(:,:,1)*pixelsize^2/(2*pi*dist*lambda);
    phi2     = phi(:,:,2)*pixelsize^2/(2*pi*dist*lambda);
    phi12    = phi1+phi2;
    maxphi1  = max(phi1(:));
    minphi1  = min(phi1(:));
    maxphi2  = max(phi2(:));
    minphi2  = min(phi2(:));
    maxphi12 = max(phi12(:));
    minphi12 = min(phi12(:));
    if errfun==1,
    errmap1  = abs((phase_object).^normpow-(phi1).^normpow);
    errmap2  = abs((phase_object).^normpow-(phi12).^normpow);
    elseif errfun==2,
    errmap1  = abs((phase_object./(pmax-pmin)).^normpow-(phi1 ./(maxphi1 -minphi1 )).^normpow);
    errmap2  = abs((phase_object./(pmax-pmin)).^normpow-(phi12./(maxphi12-minphi12)).^normpow);
    elseif errfun==3,
    errmap1  = abs(normat(phase_object).^normpow-normat(phi1).^normpow);
    errmap2  = abs(normat(phase_object).^normpow-normat(phi12).^normpow);
    end;
    err = [mean(errmap1(:)),mean(errmap2(:)),max(errmap1(:)),max(errmap2(:))];
    if ii==p0,
        intstack  = int;
        phi1stack = phi1;
        phi2stack = phi2;
        err1stack = errmap1;
        err2stack = errmap2;
        merr      = err;
        fprintf(1,['   z       [intensity]     meanLO    meanNLO (NLO/LO)      maxLO     ' ...
                   'maxNLO (NLO/LO) >phi1+phi2<     >phi1<     >phi2<\n']);
        fprintf(1,['%4.3g [%7.5g,%7.5g] %10.5g %10.5g (%5.1f%%) ' ...
                   '%10.5g %10.5g (%5.1f%%) %10.5g %10.5g %10.5g\n'], ...
                dist,min(int(:)),max(int(:)),err(1:2),100*err(2)/err(1),err(3:4),100*err(4)/err(3), ...
                maxphi12-minphi12,maxphi1-minphi1,maxphi2-minphi2);
    else,
        intstack  = cat(3,intstack,int);
        phi1stack = cat(3,phi1stack,phi1);
        phi2stack = cat(3,phi2stack,phi2);
        err1stack = cat(3,err1stack,errmap1);
        err2stack = cat(3,err2stack,errmap2);
        merr      = cat(1,merr,err);
        fprintf(1,['%4.3g [%7.5g,%7.5g] %10.5g %10.5g (%5.1f%%) ' ... 
                   '%10.5g %10.5g (%5.1f%%) %10.5g %10.5g %10.5g\n'], ...
                dist,min(int(:)),max(int(:)),err(1:2),100*err(2)/err(1),err(3:4),100*err(4)/err(3), ...
                maxphi12-minphi12,maxphi1-minphi1,maxphi2-minphi2);
    end;
end;
domain(int);