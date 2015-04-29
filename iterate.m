function iterate(dat,sli,pad,alpha);

%get uniterated phase
[phi0,phi11,phi12,phi13] = rec(dat,sli,pad,alpha);

phi = phi0 + phi11 + phi12 + phi13;

%filters:
[xi,eta]   = meshgrid(-1/2:1/(dimx):1/2-1/(dimx),-1/2:1/(dimy):1/2-1/(dimy));
xi         = fftshift(xi);
eta        = fftshift(eta);
lap        = xi.^2 + eta.^2;
inv_lap    = 1./(lap + (10^-alpha));

