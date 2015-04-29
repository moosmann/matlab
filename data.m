function [ph,rph,int,recph] = data(phtype,phdim,Ntheta,ddim);
%supply intensity pattern of a simulated phantom data. The intensity
%pattern is computed by Radon transorm of the phantom data and an
%afterward convolution of the Fresnel propagator and the transmission
%function. and so on bla bla bla


%distance of propagation between sample and screen
d      = .1; %in meter
%incident wavelength of the monochromatic beam
lambda = 0.1*10^(-9); %in meter;
%characteristic (upper) length scale of the aperture
%e.g. radius, sqrt(x*y), etc
a      = 100*10^(-6);%in meter
%Fresnel number (definition):
F      = a^2/(d*lambda);
%scaling factor
scafac = sqrt(d/lambda);
%characteristic length scale 
chalen = sqrt(lambda*d);
%detector pixel size
dx     = 1*10^(-6);%meter
dy     = dx; %meter
%detector resolution
Nx     = ddim;
Ny     = ddim;
%object function range
obj    = 10^(-9);

fprintf(1,['propagation distance: %2.3g m \n' ...
           'incident wavelength: %2.3g m \n' ... 
           'object size: %2.3g m \n' ...
           'Fresnel number: %2.3g \n'], ...
        d,lambda,a,F)

%coordinates are made dimensionless by dividing by chalen: x' =
%x/chalen. The dimensionfull phase function is given by multiplying the
%dimless phase function with scafac: phi = scafac*phi'.

%transition to normalized variables: xmax' = xmax/chalen, where for
%quadratic probes xmax = ymax. Here we have xmax = ymax = a. Thus:
%normalized variables x'' = x'/xmax' = (x/chalen)/(a/chalen) = x/a

%normalized and discretized (with respect to the pixel of the detector)
%variables: x_real = x_discrete/N_x*a_x =
%x_discrete*1/N_x*a/chalen*chalen. Thus: phi(x_real,y_real) =
%1/N*scafac*a/chalen*phi_discrete(x_discrete,y_discrete). 

  

%create 3D phantom data using the predefined phantom3d function
%phantom type
switch lower(phtype)
  case {1,'shepp-logan','sl','shepplogan'}
    ph = phantom3d('Shepp-Logan',phdim);
  case {2,'modified shepp-logan','msl','m','modifiedshepp-logan', ...
        'modifiedshepplogan'}
    ph = phantom3d('Modified Shepp-Logan',phdim);
  case {3,'square','cube','sq','cu'}
    ph = zeros(phdim,phdim,phdim);
    r  = ceil(phdim/4+1):floor(3*phdim/4);
    ph(r,r,r) = 1;
end;
ph     = -obj*ph;
[dim1,dim2,dim3] = size(ph);
fprintf(1,['dimension of phantom: %u %u %u \n' ...
           'min(phantom): %6.6g  max(phantom): %6.6g \n'], ... 
        size(ph),min(min(min(ph))),max(max(max(ph))));

%downsampling


%compute the phase function of the phantom using the  Radon transform function
theta  = 0:180/Ntheta:180*(1-1/Ntheta);
for ii = [1:dim3]
    [rph(:,ii,:),xy(ii,:,:)] = radon(ph(:,:,ii),theta);
   end;
[rdim1,rdim2,rdim3] = size(rph);
clear ii;
fprintf(1,['dimension of Radon transformed phantom: %u %u %u \n' ...
           'min(RTed phantom): %6.6g  max(RTed phantom): %6.6g \n'], ... 
        rdim1,rdim2,rdim3,min(min(min(rph))),max(max(max(rph))));
rph      = padarray(rph,[0 floor((rdim1-dim2)/2) 0],0,'pre');
rph      = padarray(rph,[0 ceil((rdim1-dim2)/2) 0],0,'post');
fprintf(1,['dimension of PADDED Radon transformed phantom: %u %u %u \n'], ...
        size(rph));


%intensity at distance d
ddim1    = ddim;
ddim2    = ddim;

if 0
    x = (ddim-rdim1)/2;
    rph = padarray(rph,[ceil(x) ceil(x) 0],0,'pre');
rph = padarray(rph,[floor(x) floor(x) 0],0,'post');
rph = fftshift(rph);
imshow(renorm(rph(:,:,1)));
end;


[xi,eta] = meshgrid(-1/2:1/(ddim2):1/2-1/(ddim2),-1/2:1/(ddim1):1/2-1/ ...
                    (ddim1));
xi       = fftshift(xi);
eta      = fftshift(eta);
fprop    = repmat(exp(-i*pi*(chalen/dx)^2*(xi.^2+eta.^2)),[1 1 length(theta)]);
ftrans   = fft2(exp(-i*2*pi*scafac*a/chalen/1*rph),ddim1,ddim2);
int      = abs(ifft2(ftrans.*fprop)).^2;

fprintf(1,['\ndimension of intensity pattern: %u %u \n' ...
           'min(intensity): %6.6g  max(intensity): %6.6g \n'], ... 
        size(int),min(min(min(int))),max(max(max(int))));

int      = int(1:rdim1,1:rdim1,:);


%3D reconstruction: Inverse RT (filtered backprojection)
if 0
    for ii = [1:ddim2]
    recph(:,:,ii) = iradon(squeeze(int(:,ii,:)),theta,'linear','Ram-Lak');
    end;
recph         = recph(2:end-1,2:end-1,:);
recphmin      = min(min(min(recph)));
recphmax      = max(max(max(recph)));
fprintf(1,['dimension of reconstructed phantom: %u %u %u \n' ...
        'min(rec. phantom): %6.6g max(rec. phantom): %6.6g \n'], ...
        size(recph),recphmin,recphmax);
recph         = (recph - recphmin)/(recphmax - recphmin);
end;


