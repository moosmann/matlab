function [m,mf_abs] = Grating(xspec,yspec,lambda,ref_ab_ind,thickness)
% Simulate grating in 1D or 2D for specified parameters: xspec = [grating
% period, sampling size], ysepc analog, thickness (of binary grating)=
% [d0,d1], ref_ab_ind = [refractive index, absorptive index], lambda =
% wavelength

    if (nargin<2) || isempty(yspec)
        yspec = xspec;
    end;
    if (nargin<3) || isempty(lambda)
        lambda = 0.496e-10;  %25 keV
    end;
    if (nargin<4) || isempty(ref_ab_ind)
        % Grating material  SU8.
        ref_ab_ind = [0,3.915e-7]; %25 keV
        %beta=6.816e-10;%25 keV
    end;
    if (nargin<5) || isempty(thickness)
        thickness = [0,lambda/ref_ab_ind(1)/2]; % pi-shift grating
    end;

    % Parameters.
    dx      = xspec(1);
    samplex = xspec(2);
    dy      = yspec(1); 
    sampley = yspec(2);
    d0      = thickness(1);
    d1      = thickness(2);
    delta   = ref_ab_ind(1);
    beta    = ref_ab_ind(2);
    [x,y]   = meshgrid(1:sampley,1:samplex);
    % Create grid.
    m = ones(samplex,sampley)*exp(2*pi*(i*delta-beta)*d0/lambda);
    m(sin(pi*(x+dx/samplex))>0) = exp(2*pi*(i*delta-beta)*d1/lambda);
    mf_abs = abs(fft(m,[],2));
    
    