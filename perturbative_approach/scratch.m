% 1,-1 grid.
pmgrid = ones(dim+1);
pmgrid(1:2:end)=-1;
pmgrid=pmgrid(1:end-1,1:end-1);
fu=zeros(dim);
n = floor(dim/4);
fu(1+n:end-n,1+n:end-n)=pmgrid(1+n:end-n,1+n:end-n);
ifft2(fu);

% LOAD PHASE MAP.
% Read Bronnikov phantom.
if 0
cd ~/data/phantom;
object = (double(mexVolRead('phase',[256 256 360],'float32')));
object = -1e2*object(:,:,1);
object = object-min(object(:));
domain(object);
u = (fft2(exp(i*object)));
ishow(log(abs(u)));
end;
% Read Shepp-Logan phantom.
cd ~/data/shepp-logan_interior_2;
[~,~,object] = loadgpbin(['sl_E24keV_z00100cm_res512_os1_000-' ...
                    'usample-real.gpbin']);
domain(object);
% Get dimensions.
[dimx,dimy]=size(object);

%display(sum(sum(exp(i*object))));
fu       = fft2(exp(1i*object),padx,pady);
% alternating grid.
pmgrid   = ones(padx+1);
pmgrid(1:2:end) = 0;
pmgrid   = pmgrid(1:end-1,1:end-1);
fph      = zeros(padx);
n        = floor(padx/4);
fph(1+n:end-n,1+n:end-n) = pmgrid(1+n:end-n,1+n:end-n);
fu       = fft2(exp(1i*ifft2(fph)),padx,pady);
% rectangle grid
fph      = zeros(padx);
n        = floor(padx/4);
fph(1+n:end-n,1+n:end-n) = 1;
fu       = fft2(exp(1i*abs(ifft2(fph))),padx,pady);

