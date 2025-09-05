p = '/data/hereon/wp/group/phase/screw';
s = dir([p '/noisy*.tif']);
nl = numel(s);

dx = 0.16e-6; %m
energy = 50000; %eV
l = E_to_lambda(energy);
z = 1.86; %m
N = 1024;
fr = dx^2 / l / z;
fprintf('\n Fresnel number : %.6f',fr)
%%

rps = -10:0.2:10;
M = numel(rps);
phast = zeros([3,numel(s),M]);
fl = cell([1,nl]);

for n = 1:nl
    %fn = [p filesep 'intensity_fr_0.0012.tif'];
    name = s(n).name;
    fn = [p filesep name];
    im = imread(fn);
    [i1] = regexp(name,'\d*');
    frn = name(i1(1):end-4);
    fl{n} = frn;
    f = str2double(frn);
    fprintf('\n name : %s',name)
    fprintf('\n noise level : %.5f',f)


    num = regexprep(sprintf('%3.1f',f),'\.','p');
    pout = [p '/pha' frn];
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
        fn = sprintf('%s/pha_noise%s_%03u_%s.tif',pout,num,m,str);
        write32bitTIFfromSingle(fn,pha)

    end    
    fprintf('\n')
end
fprintf('\n')

mopt = 62; % background flat
rp = rps(mopt);

%%
material(1).name = 'Mg10Gd';
material(1).delta =  1.29e-7;
material(1).beta = 5.88e-10;
material(2).name = 'cortical bone';
material(2).delta = 1.29e-8;
material(2).beta = 1.23e-11;

for n = 1:2
    fprintf('\n');
    m = material(n);
    m.db = m.delta/m.beta;
    m.rp = log10(m.db);
    fprintf('\n%s',m.name)
    fprintf('\n delta/beta : %.1f',m.db)
    fprintf('\n regpar : %.1f',m.rp)
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

% p = '/data/hereon/wp/group/phase/bubble/reco_tf';
% pha0 = imread('/data/hereon/wp/group/phase/bubble/sphere_phase_360.tif');
% d = dir([p '/sphere_phase*']);
% pha = read_images_to_stack(p,1,'*phase*');
% pha_mse = zeros([size(pha),numel(d)]);
% pha_nmse = zeros([size(pha),numel(d)]);
% for n = size(pha,3)
%     im = pha(:,:,n);
%     im2 = abs(im - pha0).^2;
%     pha_mse(:,:,n) = im2;
%     pha_nmse(:,:,n) = im2/mean2(im2);
% end
% 
% figure('Name','MSE')
% for n = 1:size(pha,3)
%     subplot(2,3,n)
%     imsc(pha_mse(:,:,n))
%     colorbar
% end
% 
% figure('Name','NMSE')
% for n = 1:size(pha,3)
%     subplot(2,3,n)
%     imsc(pha_nmse(:,:,n))
%     colorbar
% end