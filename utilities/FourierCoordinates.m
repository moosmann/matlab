function [xi,eta,s] = FourierCoordinates(n,ns,smax)
    
    
[xi eta s] = meshgrid(-1/2:1/n:1/2-1n,-1/2:1/n:1/2-1/n,smax/ns:smax/ns:smax);
xi       = fftshift(xi);
eta      = fftshift(eta);
k2       = xi.^2 + eta.^2;
