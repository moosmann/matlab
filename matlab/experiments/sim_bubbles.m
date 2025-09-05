%% int_invalpha10
clear all
ca

db = 2;% corresponding to db = 10^r with r = 0.3, n = 53
db = 1000;% corresponding to db = 10^r with r = 3.2, n = 103, F=0.12

N = 512;
rps = -7:0.1:7;
p = '/data/hereon/wp/group/phase/bubble/int_invalpha1000';
s = dir([p '/holo*.tiff']);
M = numel(rps);
phast = zeros([3,numel(s),M]);
fl = {'0.12','0.024','0.012','0.0024','0.0012'};

for n = 1:numel(s)
    
    name = s(n).name;
    fn = [p filesep name];
    im = imread(fn);
    fln = fl{n};
    f = str2double(fln);
    fprintf('\n name : %s',name)
    fprintf('\n F : %.5f',f)


    num = regexprep(sprintf('%7.5f',f),'\.','p');
    pout = [p '_pha' num];
    CheckAndMakePath(pout)

    imp = padarray(im,[N N],'symmetric','post');
    imf = fft2(imp);
    %phas = zeros([N,N,M],'single');
    for m = 1:M
        rp = rps(m);
        [pf, str] = PhaseFilter('tie',2*[N N],1i*f,rp);
        imp = ifft2(pf.*imf);
        pha = imp(1:N,1:N);
        phast(1,n,m) = mean2(pha);
        phast(2,n,m) = std2(pha);
        phast(3,n,m) = entropy(double(pha));
        %phas(:,:,m) = pha;
        fn = sprintf('%s/pha_%03u_%s.tif',pout,m,str);
        write32bitTIFfromSingle(fn,pha)

    end
    fprintf('\n')
end
fprintf('\n')


%% Figures
figname = 'metrics_vs_regpar';
f = figure('Name',figname);
nv = @(x) (x - min(x))./(max(x) - min(x));

subplot(1,3,1)
y = squeeze(phast(1,:,:))';
yn = nv(y);
plot(yn,'.')
title('mean')
legend(fl)

subplot(1,3,2)
y = squeeze(phast(2,:,:))';
yn = nv(y);
plot(yn,'.')
title('std')
legend(fl)


subplot(1,3,3)
y = squeeze(phast(3,:,:))';
yn = nv(y);
plot(yn,'.')
title('entropy')
legend(fl)

pfig = sprintf('%s/figures',pout);
CheckAndMakePath(pfig)
fn_fig = sprintf('%s/%s',pfig,figname);
saveas(f,fn_fig);
drawnow
%% DH data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dx = 6e-7; %m
energy = 50000; %eV
l = E_to_lambda(energy);
z = 0.0605; %m
z = 12.1;
N = 512;
fr = dx^2 / l / z;

p = '/data/hereon/wp/group/phase/bubble/';
s = dir([p '/intensity_fr_*.tif']);
rps = -10:0.2:10;
M = numel(rps);
phast = zeros([3,numel(s),M]);
fl = cell([1,4]);

for n = 1:numel(s)
    %fn = [p filesep 'intensity_fr_0.0012.tif'];
    name = s(n).name;
    fn = [p filesep name];
    im = imread(fn);
    [i1] = regexp(name,'\d*');
    frn = name(i1(1):end-4);
    fl{n} = frn;
    f = str2double(frn);
    fprintf('\n name : %s',name)
    fprintf('\n F : %.5f',f)


    num = regexprep(sprintf('%7.5f',f),'\.','p');
    pout = [p 'pha' num];
    CheckAndMakePath(pout)

    imp = padarray(im,[N N],'symmetric','post');
    imf = fft2(imp);
    %phas = zeros([N,N,M],'single');
    for m = 1:M
        rp = rps(m);
        [pf, str] = PhaseFilter('tie',2*[N N],1i*fr,rp);
        imp = ifft2(pf.*imf);
        pha = imp(1:N,1:N);
        phast(1,n,m) = mean2(pha);
        phast(2,n,m) = std2(pha);
        phast(3,n,m) = entropy(double(pha));
        %phas(:,:,m) = pha;
        fn = sprintf('%s/pha_%03u_%s.tif',pout,m,str);
        write32bitTIFfromSingle(fn,pha)

    end
    fprintf('\n')
end
fprintf('\n')


%% Figures
figname = 'metrics_vs_regpar';
f = figure('Name',figname);
nv = @(x) (x - min(x))./(max(x) - min(x));

subplot(1,3,1)
y = squeeze(phast(1,:,:))';
yn = nv(y);
plot(yn,'.')
title('mean')
legend(fl)

subplot(1,3,2)
y = squeeze(phast(2,:,:))';
yn = nv(y);
plot(yn,'.')
title('std')
legend(fl)


subplot(1,3,3)
y = squeeze(phast(3,:,:))';
yn = nv(y);
plot(yn,'.')
title('entropy')
legend(fl)

pfig = sprintf('%s/figures',pout);
CheckAndMakePath(pfig)
fn_fig = sprintf('%s/%s',pfig,figname);
saveas(f,fn_fig);
drawnow

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = '/data/hereon/wp/group/phase/bubble/reco_tf';
pha0 = imread('/data/hereon/wp/group/phase/bubble/sphere_phase_360.tif');
d = dir([p '/sphere_phase*']);
pha = read_images_to_stack(p,1,'*phase*');
pha_mse = zeros([size(pha),numel(d)]);
pha_nmse = zeros([size(pha),numel(d)]);
for n = size(pha,3)
    im = pha(:,:,n);
    im2 = abs(im - pha0).^2;
    pha_mse(:,:,n) = im2;
    pha_nmse(:,:,n) = im2/mean2(im2);
end

figure('Name','MSE')
for n = 1:size(pha,3)
    subplot(2,3,n)
    imsc(pha_mse(:,:,n))
    colorbar
end

figure('Name','NMSE')
for n = 1:size(pha,3)
    subplot(2,3,n)
    imsc(pha_nmse(:,:,n))
    colorbar
end

%%
clear all

p = '/data/hereon/wp/group/phase/bubble/sphere_phase_360.tif';

pha = imread(p);
%pha = normat(pha);
%pha = 4.7 * pha -4.4;

domain(pha)

obj = pha;
int(:,:,1) = Propagation(obj,[15 0.001 0.1e-6],2);
int(:,:,2) = Propagation(obj,[15 0.005 0.1e-6],2);
int(:,:,3) = Propagation(obj,[15 0.01 0.1e-6],2);
int(:,:,4) = Propagation(obj,[15 0.05 0.1e-6],2);
int(:,:,5) = Propagation(obj,[15 0.1 0.1e-6],2);

obj = pha - 1i*1/1000*(pha - max2(pha));
int2(:,:,1) = Propagation(obj,[15 0.001 0.1e-6],2);
int2(:,:,2) = Propagation(obj,[15 0.005 0.1e-6],2);
int2(:,:,3) = Propagation(obj,[15 0.01 0.1e-6],2);
int2(:,:,4) = Propagation(obj,[15 0.05 0.1e-6],2);
int2(:,:,5) = Propagation(obj,[15 0.1 0.1e-6],2);

%domain(int(:,:,1))

figure('Name','F = 0.12')
subplot(1,3,1)
%imsc(pha)
imagesc(pha,[0-4.4 0.3])
colorbar
axis equal tight

subplot(1,3,2)
%imsc(int(:,:,1))
im = int(:,:,1);
imagesc(im,[0.5 1.8])
colormap(gray)
colorbar
axis equal tight

subplot(1,3,3)
%imsc(int(:,:,1))
im = int2(:,:,1);
imagesc(im,[0.5 1.8])
colormap(gray)
colorbar
axis equal tight

figure('Name','F=0.024, 0.012, 0.0024, 0.0012')
for n = 1:4
    subplot(2,4,n)
    im = int(:,:,n+1);
    %imsc(im)
    %imagesc(varargin{:},'Interpolation','bilinear')
    imagesc(im,[0.5 1.8])
    colormap(gray)
    colorbar
    axis equal tight

    subplot(2,4,n+4)
    im2 = int2(:,:,n+1);
    imagesc(im,[0.5 1.8])
    colormap(gray)
    colorbar
    axis equal tight
end


