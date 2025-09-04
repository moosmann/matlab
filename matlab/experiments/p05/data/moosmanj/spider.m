clear all
ca
p = '/data/hereon/wp/user/moosmanj/spider_hair/flat_field_corr';
pout = '/data/hereon/wp/group/phase/spider_hair';
ppha = sprintf('%s/phase_maps',pout);

%proj = read_images_to_stack(p,10,'spider*tif',[],1,1);

dx = 26.2e-9; %m
energy = 11000;  %eV
l = E_to_lambda(energy);
z = 0.0789; %m
N = 2048;

fr = dx^2 / l / z;

rps = -5:0.2:7;
M = numel(rps);
phast = zeros([3,M]);

%name = 'proj_00001.tiff';
%fn = [p filesep name];
fn = '/data/hereon/wp/group/phase/spider_hair/spider_hair.tiff';
im = imread(fn);
im = rot90(im);
fprintf('\n name : %s',fn)
fprintf('\n F : %.3g',fr)
%%
CheckAndMakePath(pout)
CheckAndMakePath(ppha)

imp = padarray(im,[N N],'symmetric','post');
imf = fft2(imp);
phas = zeros([N,N,M],'single');
for m = 1:M
    rp = rps(m);
    [pf, str] = PhaseFilter('tie',2*[N N],1i*fr,rp);
    imp = ifft2(pf.*imf);
    pha = imp(1:N,1:N);
    phas(:,:,m) = pha;
    phast(1,m) = mean2(pha);
    phast(2,m) = std2(pha);
    phast(3,m) = entropy(double(pha));
    %phas(:,:,m) = pha;
    fn = sprintf('%s/pha_%02u_%s.tif',ppha,m,str);
    write32bitTIFfromSingle(fn,pha)
end
fprintf('\n')

%% Phase map
n = 33;
rp = 1.4;

name = sprintf('phase_tie_regpar%3.1f',rps(n));
name = regexprep(name,'\.','p');
fn = sprintf('%s/%s.tif',pout,name);
pha = phas(:,:,n);
write32bitTIFfromSingle(fn,pha);
fn = sprintf('%s/%s_bin2.tif',pout,name);
write32bitTIFfromSingle(fn,Binning(pha,2)/4);
fn = sprintf('%s/%s_bin4.tif',pout,name);
write32bitTIFfromSingle(fn,Binning(pha,4)/16);

%% Figures
figname = 'metrics_vs_regpar';
f = figure('Name',figname);

subplot(1,3,1)
y = squeeze(phast(1,:));
plot(y)
title('mean')

subplot(1,3,2)
y = squeeze(phast(2,:));
plot(y)
title('std')

subplot(1,3,3)
y = squeeze(phast(3,:));
plot(y)
title('entropy')

pfig = sprintf('%s/figures',pout);
CheckAndMakePath(pfig)
fn_fig = sprintf('%s/%s',pfig,figname);
saveas(f,fn_fig);
drawnow

