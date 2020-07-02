% DPC characterizationxs

%% CdWO4
font_size = 18;
line_width = 6;
prefix_list = { ...
    'ximea_cdwo4_08_028ms', ...
    'ximea_cdwo4_07_038ms', ...
    'ximea_cdwo4_06_050ms', ...
    'ximea_cdwo4_05_072ms', ...
    'ximea_cdwo4_04_102ms', ...
    'ximea_cdwo4_03_148ms', ...
    'ximea_cdwo4_02_255ms', ...
    'ximea_cdwo4_01_520ms'};
datadir = '/asap3/petra3/gpfs/p07/2019/data/11007902/raw/';
roi = [3450,3750,1300,1600]; % XIMEA % cropping roi [x1,x2,y1,y2]
vis_roi = [1,roi(2)-roi(1)+1,1,roi(4)-roi(3)+1]; % XIMEA
ims = zeros(roi(4)-roi(3)+1,roi(2)-roi(1)+1,length(prefix_list),num_img);
tic
parfor ap = 1:length(prefix_list)
    prefix = prefix_list{ap};
    imgloop = dir([datadir prefix '/*img*']);
    darkname = dir([datadir prefix '/*dar*']);
    dark = imread([datadir prefix '/' darkname.name],...
            'PixelRegion',{[roi(3) roi(4)],[roi(1) roi(2)]});
    for dpc_step = 1:num_img
        this_im = imread([datadir prefix '/' imgloop(dpc_step).name],...
            'PixelRegion',{[roi(3) roi(4)],[roi(1) roi(2)]});
        ims(:,:,ap,dpc_step) = double(this_im)-double(dark);
    end
end
fprintf([num2str(length(prefix_list)*num_img) ' images loaded, '])
toc
% get exposure time
apsize = zeros(1,length(prefix_list));
for ap = 1:length(prefix_list)
    prefix = prefix_list{ap};
    apsize(ap) = str2double(prefix(13:14));
end
exptime = zeros(1,length(prefix_list));
for et = 1:length(prefix_list)
    prefix = prefix_list{et};
    exptime(et) = str2double(prefix(16:18));
end
meancounts = zeros(1,length(prefix_list));
for ap = 1:length(prefix_list)
    these_ims = squeeze(ims(:,:,ap,:));
    tmp = these_ims(vis_roi(3):vis_roi(4),vis_roi(1):vis_roi(2),:);
    meancounts(ap) = mean(tmp(:));
end
countspers = meancounts./exptime;
% calculate visibility
vis_map = zeros(roi(4)-roi(3)+1,roi(2)-roi(1)+1,length(prefix_list));
for ap = 1:length(prefix_list)
    imsteps = squeeze(ims(:,:,ap,:));
    ft_imsteps = fft(imsteps,[],3);
    vis_map(:,:,ap) = 2*abs(ft_imsteps(:,:,2))./abs(ft_imsteps(:,:,1));
end
for ap = 1:length(prefix_list)
    vis_map(:,:,ap) = medfilt2(squeeze(vis_map(:,:,ap)),[5 5]);
end
vis_apsize = squeeze(mean(mean(vis_map(vis_roi(3):vis_roi(4),...
    vis_roi(1):vis_roi(2),:),1),2));
apsize_cdwo4 = apsize;
vis_apsize_cdwo4 = vis_apsize;
countspers_cdwo4 = countspers;
%% LuAG
datadir = '/asap3/petra3/gpfs/p07/2019/data/11007902/raw/';
prefix_list = { ...
    'ximea_luag_08_57ms',...
    'ximea_luag_07_73ms',...
    'ximea_luag_06_100ms',...
    'ximea_luag_05_140ms',...
    'ximea_luag_04_195ms',...
    'ximea_luag_03_290ms',...
    'ximea_luag_02_500ms',...
    'ximea_luag_01_1000ms'};
roi = [3450,3750,1300,1600]; % XIMEA % cropping roi [x1,x2,y1,y2]
vis_roi = [1,roi(2)-roi(1)+1,1,roi(4)-roi(3)+1]; % XIMEA
ims = zeros(roi(4)-roi(3)+1,roi(2)-roi(1)+1,length(prefix_list),num_img);
tic
parfor ap = 1:length(prefix_list)
    prefix = prefix_list{ap};
    imgloop = dir([datadir prefix '/*img*']);
    darkname = dir([datadir prefix '/*dar*']);
    dark = imread([datadir prefix '/' darkname.name],...
            'PixelRegion',{[roi(3) roi(4)],[roi(1) roi(2)]});
    for dpc_step = 1:num_img
        this_im = imread([datadir prefix '/' imgloop(dpc_step).name],...
            'PixelRegion',{[roi(3) roi(4)],[roi(1) roi(2)]});
        ims(:,:,ap,dpc_step) = double(this_im)-double(dark);
    end
end
fprintf([num2str(length(prefix_list)*num_img) ' images loaded, '])
toc
% get exposure time
apsize = zeros(1,length(prefix_list));
for ap = 1:length(prefix_list)
    prefix = prefix_list{ap};
    apsize(ap) = str2double(prefix(12:13));
end
exptime = zeros(1,length(prefix_list));
for et = 1:length(prefix_list)
    prefix = prefix_list{et};
%     exptime(et) = str2double(prefix(16:18));
    if isnan(str2double(prefix(18)))
        if isnan(str2double(prefix(17)))
            exptime(et) = str2double(prefix(15:16));
        else
            exptime(et) = str2double(prefix(15:17));
        end
    else
        exptime(et) = str2double(prefix(15:18));
    end
end
meancounts = zeros(1,length(prefix_list));
for ap = 1:length(prefix_list)
    these_ims = squeeze(ims(:,:,ap,:));
    tmp = these_ims(vis_roi(3):vis_roi(4),vis_roi(1):vis_roi(2),:);
    meancounts(ap) = mean(tmp(:));
end
countspers = meancounts./exptime;
% calculate visibility
vis_map = zeros(roi(4)-roi(3)+1,roi(2)-roi(1)+1,length(prefix_list));
for ap = 1:length(prefix_list)
    imsteps = squeeze(ims(:,:,ap,:));
    ft_imsteps = fft(imsteps,[],3);
    vis_map(:,:,ap) = 2*abs(ft_imsteps(:,:,2))./abs(ft_imsteps(:,:,1));
end
for ap = 1:length(prefix_list)
    vis_map(:,:,ap) = medfilt2(squeeze(vis_map(:,:,ap)),[5 5]);
end
vis_apsize = squeeze(mean(mean(vis_map(vis_roi(3):vis_roi(4),...
    vis_roi(1):vis_roi(2),:),1),2));
apsize_luag = apsize;
vis_apsize_luag = vis_apsize;
countspers_luag = countspers;
exptime_luag = exptime;

% plots
figure, plot(apsize_luag*0.1,vis_apsize_luag.*sqrt(countspers_luag)','.-')
hold on
plot(apsize_cdwo4*0.1,vis_apsize_cdwo4.*sqrt(countspers_cdwo4)','.-')
xlabel('aperture')
ylabel('V*\sqrt{counts/s}')
legend({'LuAG 50 \mum','CdWO_4 100 \mum'},'Location','southeast')
xlim([0,1])
title('figure of merit')

%% Visibility and Counts vs Aperture
fig = figure('Name', 'visibility | counts vs aperture', 'Units', 'Normalized', 'Position', [0.1 0.1 0.5 0.8] );
yyaxis left
p2 = plot(apsize_luag*0.1,vis_apsize_luag,'x-','Color',[0, 0.4470, 0.7410]);
hold on
p1 = plot(apsize_cdwo4*0.1,vis_apsize_cdwo4,'x-','Color',[0.8500, 0.3250, 0.0980]);
ylabel('visibility', 'FontSize',font_size )
xlabel('aperture size', 'FontSize',font_size )
yyaxis right
p3 = plot(apsize_luag*0.1,countspers_luag,'*--','Color',[0, 0.4470, 0.7410]);
hold on
p4 = plot(apsize_cdwo4*0.1,countspers_cdwo4,'*--','Color',[0.8500, 0.3250, 0.0980]);
ylabel('counts/ms', 'FontSize',font_size )
xlim([0,1])
ylim([0,90])
l1 = plot([NaN,NaN],'-','color',[0, 0.4470, 0.7410]);
l2 = plot([NaN,NaN],'-','color',[0.8500, 0.3250, 0.0980]);
legend([l1, l2],{'LuAG 50 \mum','CdWO_4 100 \mum'},'Location','southeast','FontSize',font_size )
set( p1 ,'LineWidth', line_width );
set( p2 ,'LineWidth', line_width );
set( p3 ,'LineWidth', line_width );
set( p4 ,'LineWidth', line_width );
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
ax.FontSize = font_size; 
title( fig.Name )
axis tight
filename = '/asap3/petra3/gpfs/p07/2019/data/11007902/processed/images/vis_and_counts_vs_aperture__50luag_100cdwo4.png';
saveas( fig, filename );

%% Visibilty vs Distance
datadir = '/asap3/petra3/gpfs/p07/2019/data/11007902/raw/';
prefix = 'ximea_vis12_aperture03_z';%'ximea_vis02_z';%'vis01_z'; %'ximea_vis01_z';
outerloop = dir([datadir prefix '*']);
tmp2 = dir([datadir outerloop(1).name '/*img*']);
% get sizes
zpositions = numel(outerloop); % number of z positions
num_img = numel(tmp2); % number of phase steps
% load data
refname = dir([datadir outerloop(1).name '/*ref*']);
tmp = imread([datadir outerloop(1).name '/' refname(1).name]);
%[sy,sx,sz] = size(tmp);
roi = [3750,4050,1350,1650]; % XIMEA % cropping roi [x1,x2,y1,y2]
figure, imagesc(tmp), axis equal
hold on
rectangle('Position',[roi(1),roi(3),roi(2)-roi(1),...
    roi(2)-roi(1)],'EdgeColor','r')
title('crop roi')
ims = zeros(roi(4)-roi(3)+1,roi(2)-roi(1)+1,zpositions,num_img);
tic
parfor zp = 1:zpositions
    imgloop = dir([datadir outerloop(zp).name '/*img*']);
    darkname = dir([datadir outerloop(zp).name '/*dar*']);
    dark = imread([datadir outerloop(zp).name '/' darkname.name],...
            'PixelRegion',{[roi(3) roi(4)],[roi(1) roi(2)]});
    for dpc_step = 1:num_img
        tmp = imread([datadir outerloop(zp).name '/' imgloop(dpc_step).name],...
            'PixelRegion',{[roi(3) roi(4)],[roi(1) roi(2)]});
        ims(:,:,zp,dpc_step) = double(tmp-dark);
    end
end
fprintf([num2str(zpositions*num_img) ' images loaded, '])
toc
% read distances from file names
z = zeros(1,zpositions);
for zp = 1:zpositions
    z(zp) = str2double(outerloop(zp).name(length(prefix)+1:length(prefix)+4));
end
% calculate visibility
vis_map = zeros(roi(4)-roi(3)+1,roi(2)-roi(1)+1,zpositions);
for zp = 1:zpositions
    imsteps = squeeze(ims(:,:,zp,:));
    ft_imsteps = fft(imsteps,[],3);
    vm = 2*abs(ft_imsteps(:,:,2))./abs(ft_imsteps(:,:,1));
    vis_map(:,:,zp) = medfilt2(vm,[3 3]);
end
% choose a roi to average visibility over
vis_roi = [125,175,125+25,175+25]; % XIMEA
vis_dist = squeeze(mean(mean(vis_map(vis_roi(3):vis_roi(4),...
    vis_roi(1):vis_roi(2),:),1),2));
%% Figure Distance
fig = figure('Name', 'visibility vs distance', 'Units', 'Normalized', 'Position', [0.1 0.1 0.5 0.8] );
p = plot( z, vis_dist,'.');
xlabel('z [mm]', 'FontSize',font_size )
ylabel('visibility', 'FontSize',font_size )
grid on
ylim([0,0.16])
title( fig.Name, 'FontSize',font_size )
set( p ,'LineWidth', line_width, 'MarkerSize', 14 );
ax = gca;
ax.FontSize = font_size;
legend( 'CdWO_4 100 \mum', 'Location', 'SouthEast' )
axis tight
filename = '/asap3/petra3/gpfs/p07/2019/data/11007902/processed/images/vis_vs_distance.png';
saveas( fig, filename );

%% % vis vs tilt vs distance
datadir = '/asap3/petra3/gpfs/p07/2019/data/11007902/raw/';
prefix = 'ximea_vis17_z';%'ximea_vis02_z';%'vis01_z'; %'ximea_vis01_z';
outerloop = dir([datadir prefix '*']);
z = 800:10:1000;
tiltangles = 32:1:44;
zpositions = length(z);
tilts = length(tiltangles);
ntomoscans = numel(outerloop);
num_img = numel(dir([datadir outerloop(1).name '/*img*'])); % number of phase steps
% load data
roi = [3750,4050,1000,1300]; % XIMEA % cropping roi [x1,x2,y1,y2]
vis_roi = [125,175,125+25,175+25]; % XIMEA
refname = dir([datadir outerloop(1).name '/*ref*']);
tmp = imread([datadir outerloop(1).name '/' refname(1).name]);
tmp_crop = tmp(roi(3):roi(4),roi(1):roi(2));
%[sy,sx,sz] = size(tmp);
figure, imagesc(tmp), axis equal
hold on, rectangle('Position',[roi(1),roi(3),roi(2)-roi(1),roi(2)-roi(1)],'EdgeColor','r')
title('crop roi')
figure, imagesc(tmp_crop), axis equal
hold on, rectangle('Position',[vis_roi(1),vis_roi(3),vis_roi(2)-vis_roi(1),...
    vis_roi(2)-vis_roi(1)],'EdgeColor','r')
title('vis roi')
ims = zeros(vis_roi(4)-vis_roi(3)+1,vis_roi(2)-vis_roi(1)+1,ntomoscans);
tic
parfor ii = 1:ntomoscans
    imgloop = dir([datadir outerloop(ii).name '/*img*']);
    darkname = dir([datadir outerloop(ii).name '/*dar*']);
    dark = imread([datadir outerloop(ii).name '/' darkname.name],...
            'PixelRegion',{[roi(3) roi(4)],[roi(1) roi(2)]});
    for dpc_step = 1:num_img
        tmp = imread([datadir outerloop(ii).name '/' imgloop(dpc_step).name],...
            'PixelRegion',{[roi(3) roi(4)],[roi(1) roi(2)]});
        tmp = double(tmp-dark);
        ims(:,:,ii,dpc_step) = tmp(vis_roi(3):vis_roi(4),vis_roi(1):vis_roi(2));
    end
end
fprintf([num2str(ntomoscans*num_img) ' images loaded, '])
toc
vis_map = zeros(vis_roi(4)-vis_roi(3)+1,vis_roi(2)-vis_roi(1)+1,ntomoscans);
for ii = 1:ntomoscans
    these_ims = squeeze(ims(:,:,ii,:));
    ft_ims = fft(these_ims,[],3);
    vm = 2*abs(ft_ims(:,:,2))./abs(ft_ims(:,:,1));
    vis_map(:,:,ii) = medfilt2(vm,[3 3]);
end
%%
fun = @(block_struct)  ...
    median(block_struct.data )*ones(size(block_struct.data));
vis_dat = zeros(1,ntomoscans);
for ii = 1:ntomoscans
    vm = vis_map(:,:,ii) ;
    vm2 = blockproc(vm,[5 5],fun);
    vis_dat(ii) = median(vm2,'all');
end
vis_dat2 = reshape(vis_dat,[tilts,zpositions]);
figure, plot(z,vis_dat2(1:3:end,:),'.-')
xlabel('z [mm]')
ylabel('visibility')
legend(num2str(tiltangles(1:3:end)'))
title('Visibility, distance, and tilt angle')

%% phase projections
datadir = '/asap3/petra3/gpfs/p07/2019/data/11007454/raw/';
sdir = 'bmc10_mouse21_xsgi_a/tiff0000'; % scan directory
fnames = dir([datadir sdir filesep '*img*']);
flatnames = dir([datadir sdir filesep '*ref*']);
darkname = dir([datadir sdir filesep '*dar*']);
crop_roi = [1,7920,1,2600];
dark = imread([datadir sdir filesep darkname(1).name],...
            'PixelRegion',{[crop_roi(3) crop_roi(4)],[crop_roi(1) crop_roi(2)]});
num_img = 5;
sdir = 'bmc10_mouse21_xsgi_a/tiff0010'; % scan directory
fnames = dir([datadir sdir filesep '*img*']);
flatnames = dir([datadir sdir filesep '*ref*']);
ims = zeros(crop_roi(4)-crop_roi(3)+1,crop_roi(2)-crop_roi(1)+1,num_img);
flat_ims = zeros(crop_roi(4)-crop_roi(3)+1,crop_roi(2)-crop_roi(1)+1,num_img);
tic
for ii = 1:num_img
    tmp = imread([datadir sdir filesep fnames(ii+5*500).name],...
            'PixelRegion',{[crop_roi(3) crop_roi(4)],[crop_roi(1) crop_roi(2)]});
    ims(:,:,ii) = tmp-dark;
    tmp = imread([datadir sdir filesep flatnames(ii).name],...
            'PixelRegion',{[crop_roi(3) crop_roi(4)],[crop_roi(1) crop_roi(2)]});
    flat_ims(:,:,ii) = tmp-dark;
end
toc
% dpc processing
wrap = @(datin) ((datin+pi)-2*pi*floor((datin+pi)/(2*pi)))-pi;
ft_ims = fft(ims,[],3);
ft_flat_ims = fft(flat_ims,[],3);
dp_im = wrap(angle(ft_ims(:,:,2))-angle(ft_flat_ims(:,:,2)));
df_im = (abs(ft_ims(:,:,2))./abs(ft_ims(:,:,1)))./(abs(ft_flat_ims(:,:,2))./abs(ft_flat_ims(:,:,1)));
ac_im = abs(ft_ims(:,:,1))./abs(ft_flat_ims(:,:,1));
vis_ff = 2*abs(ft_flat_ims(:,:,2))./abs(ft_flat_ims(:,:,1));

figure, imagesc(ims(:,:,1))
title('raw proj')

figure, subplot(311)
imagesc(dp_im,[-pi,pi]/5), axis equal
ylim([crop_roi(3),crop_roi(4)])
xlim([crop_roi(1),crop_roi(2)])
title('differential phase')
subplot(312), imagesc(ac_im,[0.7,1.1]), axis equal
ylim([crop_roi(3),crop_roi(4)])
xlim([crop_roi(1),crop_roi(2)])
title('attenuation')
subplot(313), imagesc(df_im,[0,2.5]), axis equal
ylim([crop_roi(3),crop_roi(4)])
xlim([crop_roi(1),crop_roi(2)])
title('dark field')

% % image of the gratings
tmp = flat_ims(:,:,1);
figure, imagesc(tmp), axis equal
ylim([crop_roi(3),crop_roi(4)])
xlim([crop_roi(1),crop_roi(2)])
title('fringe pattern at d = 900 mm')

roi = [3750,4050,1000,1300]; % XIMEA % cropping roi [x1,x2,y1,y2]
vis_roi = [125,175,125+25,175+25]; % XIMEA

figure, subplot(131)
imagesc(tmp), axis equal
ylim([crop_roi(3),crop_roi(4)])
xlim([crop_roi(1),crop_roi(2)])
hold on, rectangle('Position',[roi(1),roi(3),roi(2)-roi(1),roi(2)-roi(1)],'EdgeColor','r')
title('fringe pattern at d = 900 mm')
tmp_crop = tmp(roi(3):roi(4),roi(1):roi(2));
subplot(132), imagesc(tmp_crop), axis equal
ylim([1,300])
xlim([1,300])
hold on, rectangle('Position',[vis_roi(1),vis_roi(3),vis_roi(2)-vis_roi(1),...
    vis_roi(2)-vis_roi(1)],'EdgeColor','r')
vis_tmp_crop = tmp_crop(vis_roi(3):vis_roi(4),vis_roi(1):vis_roi(2));
subplot(133), imagesc(vis_tmp_crop), axis equal
ylim([1,50])
xlim([1,50])

%% dpc scripts
% data is in '/asap3/petra3/gpfs/p07/2019/data/11007902/raw/'
datadir = '/asap3/petra3/gpfs/p07/2019/data/11007454/raw/';
sdir1 = '/bmc13_sheep_brain_formalin_a'; % scan directory

dirs = dir([datadir sdir1 filesep 'tiff*']);
ndirs = numel(dirs);
n_proj = zeros(1,ndirs);
n_frames= zeros(1,ndirs);
for ii =1:ndirs
    sdir2 = [sdir1 filesep 'tiff' num2str(ii-1,'%04d')];
    fnames = dir([datadir sdir2 filesep '*img*']);
    flatnames = dir([datadir sdir2 filesep '*ref*']);
    darkname = dir([datadir sdir2 filesep '*dar*']);
    n_proj(ii) = numel(fnames);
    n_frames(ii) = numel(fnames)+numel(flatnames)+numel(darkname);
end

a=imread([datadir sdir2 filesep fnames(end-100).name]);
figure,imagesc(a)
%% check that scan is running
datadir = '/asap3/petra3/gpfs/p07/2019/data/11007454/raw/';
sdir = 'bmc10_mouse21_xsgi_a/tiff0000'; % scan directory
fnames = dir([datadir sdir filesep '*img*']);
flatnames = dir([datadir sdir filesep '*ref*']);
darkname = dir([datadir sdir filesep '*dar*']);
dark = single(imread([datadir sdir filesep darkname(1).name]));
im = single(imread([datadir sdir filesep fnames(end).name]))-dark;
ref = single(imread([datadir sdir filesep flatnames(end).name]))-dark;

figure, imagesc(im./ref,[0.25,1.1])

%% one phase projection
datadir = '/asap3/petra3/gpfs/p07/2019/data/11007454/raw/';
sdir = 'vis_0512_01/tiff0000'; % scan directory

fnames = dir([datadir sdir filesep '*img*']);
flatnames = dir([datadir sdir filesep '*ref*']);
darkname = dir([datadir sdir filesep '*dar*']);

crop_roi = [1,7920,1,2600];
dark = imread([datadir sdir filesep darkname(1).name],...
            'PixelRegion',{[crop_roi(3) crop_roi(4)],[crop_roi(1) crop_roi(2)]});

num_img = 5;
ims = zeros(crop_roi(4)-crop_roi(3)+1,crop_roi(2)-crop_roi(1)+1,num_img);
flat_ims = zeros(crop_roi(4)-crop_roi(3)+1,crop_roi(2)-crop_roi(1)+1,num_img);
tic
for ii = 1:num_img
    tmp = imread([datadir sdir filesep fnames(ii).name],...
            'PixelRegion',{[crop_roi(3) crop_roi(4)],[crop_roi(1) crop_roi(2)]});
    ims(:,:,ii) = tmp-dark;
    tmp = imread([datadir sdir filesep flatnames(ii).name],...
            'PixelRegion',{[crop_roi(3) crop_roi(4)],[crop_roi(1) crop_roi(2)]});
    flat_ims(:,:,ii) = tmp-dark;
end
toc

% dpc processing
wrap = @(datin) ((datin+pi)-2*pi*floor((datin+pi)/(2*pi)))-pi;

ft_ims = fft(ims,[],3);
ft_flat_ims = fft(flat_ims,[],3);

dp_im = wrap(angle(ft_ims(:,:,2))-angle(ft_flat_ims(:,:,2)));
df_im = (abs(ft_ims(:,:,2))./abs(ft_ims(:,:,1)))./(abs(ft_flat_ims(:,:,2))./abs(ft_flat_ims(:,:,1)));
ac_im = abs(ft_ims(:,:,1))./abs(ft_flat_ims(:,:,1));
vis_ff = 2*abs(ft_flat_ims(:,:,2))./abs(ft_flat_ims(:,:,1));

figure, imagesc(ims(:,:,1))
title('raw proj')

figure, imagesc(dp_im,[-pi,pi]/5), axis equal
title('dpc')

figure, imagesc(ac_im,[0.7,1.1]), axis equal
title('att')

figure, imagesc(df_im,[0,2.5]), axis equal
title('df')

figure, imagesc(ims(200:450,3200:3450,1))

figure, imagesc(dp_im(200:450,3200:3450))
figure, imagesc(wrap(angle(ft_ims(200:450,3200:3450,2))))
figure, imagesc(wrap(angle(ft_flat_ims(200:450,3200:3450,2))))

roi = [2500,2750,1100,1350];
figure, imagesc(vis_ff(roi(3):roi(4),roi(1):roi(2)))

median(vis_ff(roi(3):roi(4),roi(1):roi(2)),'all')

figure, imagesc(vis_ff,[0.07,0.17])

%% dpc reco
datadir = '/asap3/petra3/gpfs/p07/2019/data/11007902/raw/';
sdir = 'bmc01_greek_brain_e'; % scan directory

fnames = dir([datadir sdir filesep '*img*']);
flatnames = dir([datadir sdir filesep '*ref*']);
darkname = dir([datadir sdir filesep '*dar*']);

% define a cropping region
tmp = imread([datadir sdir filesep fnames(1).name]);
croi = [360,7561,1201,1401]; % [x1,x2,y1,y2]
figure, imagesc(tmp), hold on
rectangle('Position',[croi(1),croi(3),croi(2)-croi(1),...
    croi(4)-croi(3)],'EdgeColor','r')
title('crop roi')

% load dark
tmp = imread([datadir sdir filesep fnames(ii).name],...
            'PixelRegion',{[crop_roi(3) crop_roi(4)],[crop_roi(1) crop_roi(2)]});


%% one phase projection
datadir = '/asap3/petra3/gpfs/p07/2019/data/11007902/raw/';
sdir = 'dpc_proj01'; % scan directory

fnames = dir([datadir sdir filesep '*img*']);
flatnames = dir([datadir sdir filesep '*ref*']);
darkname = dir([datadir sdir filesep '*dar*']);

crop_roi = [101,7000,1,650];
dark = imread([datadir sdir filesep darkname(1).name],...
            'PixelRegion',{[crop_roi(3) crop_roi(4)],[crop_roi(1) crop_roi(2)]});

num_img = numel(fnames);
ims = zeros(crop_roi(4)-crop_roi(3)+1,crop_roi(2)-crop_roi(1)+1,num_img);
flat_ims = zeros(crop_roi(4)-crop_roi(3)+1,crop_roi(2)-crop_roi(1)+1,num_img);
tic
parfor ii = 1:num_img
    tmp = imread([datadir sdir filesep fnames(ii).name],...
            'PixelRegion',{[crop_roi(3) crop_roi(4)],[crop_roi(1) crop_roi(2)]});
    ims(:,:,ii) = tmp-dark;
    tmp = imread([datadir sdir filesep flatnames(ii).name],...
            'PixelRegion',{[crop_roi(3) crop_roi(4)],[crop_roi(1) crop_roi(2)]});
    flat_ims(:,:,ii) = tmp-dark;
end
toc

% dpc processing
wrap = @(datin) ((datin+pi)-2*pi*floor((datin+pi)/(2*pi)))-pi;

ft_ims = fft(ims,[],3);
ft_flat_ims = fft(flat_ims,[],3);

dp_im = wrap(angle(ft_ims(:,:,2))-angle(ft_flat_ims(:,:,2)));
df_im = (abs(ft_ims(:,:,2))./abs(ft_ims(:,:,1)))./(abs(ft_flat_ims(:,:,2))./abs(ft_flat_ims(:,:,1)));
ac_im = abs(ft_ims(:,:,1))./abs(ft_flat_ims(:,:,1));
vis_ff = 2*abs(ft_flat_ims(:,:,2))./abs(ft_flat_ims(:,:,1));

figure, imagesc(ims(:,:,1))
title('raw proj')

figure, imagesc(dp_im,[-pi,pi]/5), axis equal
title('dpc')

figure, imagesc(ac_im,[0.7,1.1]), axis equal
title('att')

figure, imagesc(df_im,[0,2.5]), axis equal
title('df')

figure, imagesc(ims(200:450,3200:3450,1))

figure, imagesc(dp_im(200:450,3200:3450))
figure, imagesc(wrap(angle(ft_ims(200:450,3200:3450,2))))
figure, imagesc(wrap(angle(ft_flat_ims(200:450,3200:3450,2))))

%% checking aperture size (second  scintillator)
% prefix_list = {'ximea_luag_08_57ms','ximea_luag_07_73ms',...
%     'ximea_luag_06_100ms','ximea_luag_05_140ms',...
%     'ximea_luag_04_195ms','ximea_luag_03_290ms',...
%     'ximea_luag_02_500ms','ximea_luag_01_1000ms'};
prefix_list = {'ximea_cdwo4_08_028ms','ximea_cdwo4_07_038ms',...
    'ximea_cdwo4_06_050ms','ximea_cdwo4_05_072ms',...
    'ximea_cdwo4_04_102ms','ximea_cdwo4_03_148ms',...
    'ximea_cdwo4_02_255ms','ximea_cdwo4_01_520ms'};

datadir = '/asap3/petra3/gpfs/p07/2019/data/11007902/raw/';

roi = [3450,3750,1300,1600]; % XIMEA % cropping roi [x1,x2,y1,y2]
%vis_roi = [100,150,100,150]; % XIMEA
vis_roi = [1,roi(2)-roi(1)+1,1,roi(4)-roi(3)+1]; % XIMEA

ims = zeros(roi(4)-roi(3)+1,roi(2)-roi(1)+1,length(prefix_list),num_img);
tic
parfor ap = 1:length(prefix_list)
    prefix = prefix_list{ap};
    imgloop = dir([datadir prefix '/*img*']);
    darkname = dir([datadir prefix '/*dar*']);
    dark = imread([datadir prefix '/' darkname.name],...
            'PixelRegion',{[roi(3) roi(4)],[roi(1) roi(2)]});
    for dpc_step = 1:num_img
        this_im = imread([datadir prefix '/' imgloop(dpc_step).name],...
            'PixelRegion',{[roi(3) roi(4)],[roi(1) roi(2)]});
        ims(:,:,ap,dpc_step) = double(this_im)-double(dark);
    end
end
fprintf([num2str(length(prefix_list)*num_img) ' images loaded, '])
toc

% get exposure time
apsize = zeros(1,length(prefix_list));
for ap = 1:length(prefix_list)
    prefix = prefix_list{ap};
    apsize(ap) = str2double(prefix(13:14));
end

exptime = zeros(1,length(prefix_list));
for et = 1:length(prefix_list)
    prefix = prefix_list{et};
    exptime(et) = str2double(prefix(16:18));
%     if isnan(str2double(prefix(18)))
%         if isnan(str2double(prefix(17)))
%             exptime(et) = str2double(prefix(15:16));
%         else
%             exptime(et) = str2double(prefix(15:17));
%         end
%     else
%         exptime(et) = str2double(prefix(15:18));
%     end
end

meancounts = zeros(1,length(prefix_list));
for ap = 1:length(prefix_list)
    these_ims = squeeze(ims(:,:,ap,:));
    tmp = these_ims(vis_roi(3):vis_roi(4),vis_roi(1):vis_roi(2),:);
    meancounts(ap) = mean(tmp(:));
end

countspers = meancounts./exptime;


% calculate visibility
vis_map = zeros(roi(4)-roi(3)+1,roi(2)-roi(1)+1,length(prefix_list));
for ap = 1:length(prefix_list)
    imsteps = squeeze(ims(:,:,ap,:));
    ft_imsteps = fft(imsteps,[],3);
    vis_map(:,:,ap) = 2*abs(ft_imsteps(:,:,2))./abs(ft_imsteps(:,:,1));
end

% fun = @(block_struct)  ...
%     median(block_struct.data,'all')*ones(size(block_struct.data));
% vis_dat = zeros(1,ntomoscans);
% for ii = 1:ntomoscans
%     vm = vis_map(:,:,ii) ;
%     vm2 = blockproc(vm,[5 5],fun);
%     vis_dat(ii) = median(vm2,'all');
% end

for ap = 1:length(prefix_list)
    vis_map(:,:,ap) = medfilt2(squeeze(vis_map(:,:,ap)),[5 5]);
end
vis_apsize = squeeze(mean(mean(vis_map(vis_roi(3):vis_roi(4),...
    vis_roi(1):vis_roi(2),:),1),2));

figure, plot(apsize*0.1,vis_apsize,'.-')
xlabel('aperture size')
ylabel('visibility')
xlim([0,1])
title('CdWO_4')

figure, plot(apsize*0.1,vis_apsize.*sqrt(countspers)','.-')
xlabel('aperture')
ylabel('V*sqrt(counts/s)')
title('CdWO_4')

figure, yyaxis left
plot(apsize*0.1,vis_apsize,'.-')
ylabel('visibility')
xlabel('aperture size')
yyaxis right
plot(apsize*0.1,countspers,'.-')
ylabel('counts/ms')
xlim([0,1])
title('CdWO_4')


%% compare before and after focusing
datadir = '/asap3/petra3/gpfs/p07/2019/data/11007902/raw/';
sdir = 'ximea_vis18_z0900mm_g1rotx044'; % scan directory

fnames = dir([datadir sdir filesep '*img*']);
darkname = dir([datadir sdir filesep '*dar*']);
dark = imread([datadir sdir filesep darkname(1).name]);
[sy,sx] = size(dark);

num_img = numel(fnames);
ims = zeros(sy,sx,num_img);
tic
parfor ii = 1:num_img
    tmp = imread([datadir sdir filesep fnames(ii).name]);
    ims(:,:,ii) = tmp-dark;
end
toc

figure, imagesc(ims(:,:,1))

ft_ims = fft(ims,[],3);
vis_map = abs(ft_ims(:,:,2))./abs(ft_ims(:,:,1));

figure, imagesc(vis_map,[0,12]*1e-2), colorbar
title('focused')


datadir = '/asap3/petra3/gpfs/p07/2019/data/11007902/raw/';
sdir = 'ximea_vis17_z0900mm_g1rotx044'; % scan directory

fnames = dir([datadir sdir filesep '*img*']);
darkname = dir([datadir sdir filesep '*dar*']);
dark = imread([datadir sdir filesep darkname(1).name]);
[sy,sx] = size(dark);

num_img = numel(fnames);
ims = zeros(sy,sx,num_img);
tic
parfor ii = 1:num_img
    tmp = imread([datadir sdir filesep fnames(ii).name]);
    ims(:,:,ii) = tmp-dark;
end
toc

figure, imagesc(ims(:,:,1))

ft_ims = fft(ims,[],3);
vis_map2 = abs(ft_ims(:,:,2))./abs(ft_ims(:,:,1));

figure, imagesc(vis_map2,[0,12]*1e-2), colorbar
title('before focusing')


%% 

%% visibility vs. detector position
% % Note: if you run this during acquisition, index errors just mean not
% % all projections are taken so try running it again
datadir = '/asap3/petra3/gpfs/p07/2019/data/11007902/raw/';
prefix = 'ximea_vis12_aperture03_z';%'ximea_vis02_z';%'vis01_z'; %'ximea_vis01_z';
outerloop = dir([datadir prefix '*']);
tmp2 = dir([datadir outerloop(1).name '/*img*']);

% get sizes
zpositions = numel(outerloop); % number of z positions
num_img = numel(tmp2); % number of phase steps

% load data
refname = dir([datadir outerloop(1).name '/*ref*']);
tmp = imread([datadir outerloop(1).name '/' refname(1).name]);
[sy,sx,sz] = size(tmp);
%figure, imagesc(tmp), axis equal
roi = [3750,4050,1350,1650]; % XIMEA % cropping roi [x1,x2,y1,y2]
%roi = [1750,2050,2050,2350]; % KIT % cropping roi [x1,x2,y1,y2]
figure, imagesc(tmp), axis equal
hold on
rectangle('Position',[roi(1),roi(3),roi(2)-roi(1),...
    roi(2)-roi(1)],'EdgeColor','r')
title('crop roi')

ims = zeros(roi(4)-roi(3)+1,roi(2)-roi(1)+1,zpositions,num_img);
tic
parfor zp = 1:zpositions
    imgloop = dir([datadir outerloop(zp).name '/*img*']);
    darkname = dir([datadir outerloop(zp).name '/*dar*']);
    dark = imread([datadir outerloop(zp).name '/' darkname.name],...
            'PixelRegion',{[roi(3) roi(4)],[roi(1) roi(2)]});
    for dpc_step = 1:num_img
        tmp = imread([datadir outerloop(zp).name '/' imgloop(dpc_step).name],...
            'PixelRegion',{[roi(3) roi(4)],[roi(1) roi(2)]});
        ims(:,:,zp,dpc_step) = double(tmp-dark);
    end
end
fprintf([num2str(zpositions*num_img) ' images loaded, '])
toc

% read distances from file names
z = zeros(1,zpositions);
for zp = 1:zpositions
    z(zp) = str2double(outerloop(zp).name(length(prefix)+1:length(prefix)+4));
end

% calculate visibility
vis_map = zeros(roi(4)-roi(3)+1,roi(2)-roi(1)+1,zpositions);
for zp = 1:zpositions
    imsteps = squeeze(ims(:,:,zp,:));
    ft_imsteps = fft(imsteps,[],3);
    vm = 2*abs(ft_imsteps(:,:,2))./abs(ft_imsteps(:,:,1));
    vis_map(:,:,zp) = medfilt2(vm,[3 3]);
end

% choose a roi to average visibility over
%vis_roi = [125,175,125,175]; % KIT
vis_roi = [125,175,125+25,175+25]; % XIMEA
% figure, imagesc(imsteps(:,:,1)), axis equal
% hold on
% rectangle('Position',[vis_roi(1),vis_roi(3),vis_roi(2)-vis_roi(1),...
%     vis_roi(2)-vis_roi(1)],'EdgeColor','r')

vis_dist = squeeze(mean(mean(vis_map(vis_roi(3):vis_roi(4),...
    vis_roi(1):vis_roi(2),:),1),2));

figure, plot(z,vis_dist,'.-')
xlabel('z [mm]')
ylabel('visibility')
grid on
title(prefix,'Interpreter','none')

% vis_dist_kit = vis_dist;
% z_kit = z;
% vis_dist_ximea = vis_dist;
% z_ximea = z;
% vis_dist_ximea2 = vis_dist;
% z_ximea2 = z;
% vis_dist_ximea3 = vis_dist;
% z_ximea3 = z;

figure, plot(z_kit,vis_dist_kit,'.-')
hold on
plot(z_ximea,vis_dist_ximea,'.-')
xlabel('z [mm]')
ylabel('visibility')
title('vis vs. distance')
grid on
legend({'kit','ximea'},'Location','northwest')

%% visibility vs. detector position and tilt angle
% % Note: if you run this during acquisition, index errors just mean not
% % all projections are taken so try running it again
datadir = '/asap3/petra3/gpfs/p07/2019/data/11007902/raw/';
prefix = 'ximea_vis17_z';%'ximea_vis02_z';%'vis01_z'; %'ximea_vis01_z';
outerloop = dir([datadir prefix '*']);

z = 800:10:1000;
tiltangles = 32:1:44;

zpositions = length(z);
tilts = length(tiltangles);
ntomoscans = numel(outerloop);
num_img = numel(dir([datadir outerloop(1).name '/*img*'])); % number of phase steps

% load data
roi = [3750,4050,1000,1300]; % XIMEA % cropping roi [x1,x2,y1,y2]
vis_roi = [125,175,125+25,175+25]; % XIMEA
refname = dir([datadir outerloop(1).name '/*ref*']);
tmp = imread([datadir outerloop(1).name '/' refname(1).name]);
tmp_crop = tmp(roi(3):roi(4),roi(1):roi(2));
[sy,sx,sz] = size(tmp);
figure, imagesc(tmp), axis equal
hold on, rectangle('Position',[roi(1),roi(3),roi(2)-roi(1),roi(2)-roi(1)],'EdgeColor','r')
title('crop roi')
figure, imagesc(tmp_crop), axis equal
hold on, rectangle('Position',[vis_roi(1),vis_roi(3),vis_roi(2)-vis_roi(1),...
    vis_roi(2)-vis_roi(1)],'EdgeColor','r')
title('vis roi')

ims = zeros(vis_roi(4)-vis_roi(3)+1,vis_roi(2)-vis_roi(1)+1,ntomoscans);
tic
parfor ii = 1:ntomoscans
    imgloop = dir([datadir outerloop(ii).name '/*img*']);
    darkname = dir([datadir outerloop(ii).name '/*dar*']);
    dark = imread([datadir outerloop(ii).name '/' darkname.name],...
            'PixelRegion',{[roi(3) roi(4)],[roi(1) roi(2)]});
    for dpc_step = 1:num_img
        tmp = imread([datadir outerloop(ii).name '/' imgloop(dpc_step).name],...
            'PixelRegion',{[roi(3) roi(4)],[roi(1) roi(2)]});
        tmp = double(tmp-dark);
        ims(:,:,ii,dpc_step) = tmp(vis_roi(3):vis_roi(4),vis_roi(1):vis_roi(2));
    end
end
fprintf([num2str(ntomoscans*num_img) ' images loaded, '])
toc

vis_map = zeros(vis_roi(4)-vis_roi(3)+1,vis_roi(2)-vis_roi(1)+1,ntomoscans);
for ii = 1:ntomoscans
    these_ims = squeeze(ims(:,:,ii,:));
    ft_ims = fft(these_ims,[],3);
    vm = 2*abs(ft_ims(:,:,2))./abs(ft_ims(:,:,1));
    vis_map(:,:,ii) = medfilt2(vm,[3 3]);
end

fun = @(block_struct)  ...
    median(block_struct.data,'all')*ones(size(block_struct.data));
vis_dat = zeros(1,ntomoscans);
for ii = 1:ntomoscans
    vm = vis_map(:,:,ii) ;
    vm2 = blockproc(vm,[5 5],fun);
    vis_dat(ii) = median(vm2,'all');
end

figure, plot(tiltangles,vis_dat(1:tilts),'.-')
xlabel('tilt angle [deg]')
ylabel('visibility')

vis_dat2 = reshape(vis_dat,[tilts,zpositions]);

figure, plot(tiltangles,vis_dat2,'.-')
xlabel('tilt angle [deg]')
ylabel('visibility')

figure, plot(z,vis_dat2(13,:),'.-')
hold on
plot(z,vis_dat2(7,:),'.-')
xlabel('z [mm]')
ylabel('visibility')
legend({[num2str(tiltangles(13))],[num2str(tiltangles(7))]})

figure,plot(z_ximea,vis_dist_ximea,'.-')
hold on
plot(z,vis_dat2(13,:),'.-')
plot(z,vis_dat2(7,:),'.-')
xlabel('z [mm]')
ylabel('visibility')
legend({'ximea, 38',[num2str(tiltangles(13))],[num2str(tiltangles(7))]})

figure, imagesc(z,tiltangles,vis_dat2,[0.125,0.145])
colorbar
xlabel('z [mm]')
ylabel('tilt angle [deg]')
title(prefix,'Interpreter','none')

%% visibility vs. detector position and tilt angle (LuAG)
% % Note: if you run this during acquisition, index errors just mean not
% % all projections are taken so try running it again
datadir = '/asap3/petra3/gpfs/p07/2019/data/11007902/raw/';
prefix = 'ximea_vis20_z';%'ximea_vis02_z';%'vis01_z'; %'ximea_vis01_z';
outerloop = dir([datadir prefix '*']);

z = 700:25:1100;
tiltangles = 30:2:50;

zpositions = length(z);
tilts = length(tiltangles);
ntomoscans = numel(outerloop);
num_img = numel(dir([datadir outerloop(1).name '/*img*'])); % number of phase steps

% load data
roi = [3750,4050,1000,1300]; % XIMEA % cropping roi [x1,x2,y1,y2]
vis_roi = [125,175,125+25,175+25]; % XIMEA
refname = dir([datadir outerloop(1).name '/*ref*']);
tmp = imread([datadir outerloop(1).name '/' refname(1).name]);
tmp_crop = tmp(roi(3):roi(4),roi(1):roi(2));
[sy,sx,sz] = size(tmp);
figure, imagesc(tmp), axis equal
hold on, rectangle('Position',[roi(1),roi(3),roi(2)-roi(1),roi(2)-roi(1)],'EdgeColor','r')
title('crop roi')
figure, imagesc(tmp_crop), axis equal
hold on, rectangle('Position',[vis_roi(1),vis_roi(3),vis_roi(2)-vis_roi(1),...
    vis_roi(2)-vis_roi(1)],'EdgeColor','r')
title('vis roi')

ims = zeros(vis_roi(4)-vis_roi(3)+1,vis_roi(2)-vis_roi(1)+1,ntomoscans);
tic
parfor ii = 1:ntomoscans
    imgloop = dir([datadir outerloop(ii).name '/*img*']);
    darkname = dir([datadir outerloop(ii).name '/*dar*']);
    dark = imread([datadir outerloop(ii).name '/' darkname.name],...
            'PixelRegion',{[roi(3) roi(4)],[roi(1) roi(2)]});
    for dpc_step = 1:num_img
        tmp = imread([datadir outerloop(ii).name '/' imgloop(dpc_step).name],...
            'PixelRegion',{[roi(3) roi(4)],[roi(1) roi(2)]});
        tmp = double(tmp-dark);
        ims(:,:,ii,dpc_step) = tmp(vis_roi(3):vis_roi(4),vis_roi(1):vis_roi(2));
    end
end
fprintf([num2str(ntomoscans*num_img) ' images loaded, '])
toc

vis_map = zeros(vis_roi(4)-vis_roi(3)+1,vis_roi(2)-vis_roi(1)+1,ntomoscans);
for ii = 1:ntomoscans
    these_ims = squeeze(ims(:,:,ii,:));
    ft_ims = fft(these_ims,[],3);
    vm = 2*abs(ft_ims(:,:,2))./abs(ft_ims(:,:,1));
    vis_map(:,:,ii) = medfilt2(vm,[3 3]);
end

fun = @(block_struct)  ...
    median(block_struct.data,'all')*ones(size(block_struct.data));
vis_dat = zeros(1,ntomoscans);
for ii = 1:ntomoscans
    vm = vis_map(:,:,ii) ;
    vm2 = blockproc(vm,[5 5],fun);
    vis_dat(ii) = median(vm2,'all');
end

vis_dat2 = reshape(vis_dat,[tilts,zpositions]);

figure, plot(tiltangles,vis_dat2,'.-')
xlabel('tilt angle [deg]')
ylabel('visibility')


figure,plot(z_ximea,vis_dist_ximea,'.-')
hold on
plot(z,vis_dat2,'.-')
xlabel('z [mm]')
ylabel('visibility')

figure, imagesc(z,tiltangles,vis_dat2,[0.085,0.105])
colorbar
xlabel('z [mm]')
ylabel('tilt angle [deg]')
title(prefix,'Interpreter','none')


%% checking settle times
% 22 ms exposure time
datadir = '/asap3/petra3/gpfs/p07/2019/data/11007902/raw/';
prefix_list = {'ximea_vis04_settime0000ms','ximea_vis04_settime0010ms',...
    'ximea_vis04_settime0020ms','ximea_vis04_settime0030ms','ximea_vis04_settime0040ms',...
    'ximea_vis04_settime0050ms','ximea_vis04_settime0070ms',...
    'ximea_vis04_settime0080ms','ximea_vis04_settime0090ms','ximea_vis04_settime0100ms',...
    'ximea_vis04_settime0150ms','ximea_vis04_settime0200ms',...
    'ximea_vis06_settime0020ms_zero'};
roi = [3750,4050,1350+900,1650+900]; % XIMEA % cropping roi [x1,x2,y1,y2]
vis_roi = [125,175,125+25,175+25]; % XIMEA
num_sett = length(prefix_list);
num_img = length( dir([datadir prefix '/*img*']) );
ims = zeros(roi(4)-roi(3)+1,roi(2)-roi(1)+1,num_sett,num_img);
ni = 3;
ims2 = zeros(roi(4)-roi(3)+1,roi(2)-roi(1)+1,num_sett,ni*num_img);
tic
for st = 1:num_sett
    prefix = prefix_list{st};
    fprintf( '\n %s:', prefix )
    imgloop = dir([datadir prefix '/*img*']);
    refloop = dir([datadir prefix '/*ref*']);
    darkname = dir([datadir prefix '/*dar*']);
    dark = imread([datadir prefix '/' darkname.name],...
            'PixelRegion',{[roi(3) roi(4)],[roi(1) roi(2)]});
    dark = FilterPixel( dark );
    for dpc_step = 1:num_img
        fprintf( ' %u', dpc_step )
        
        % Compare img and ref at the same position, img and ref series were
        % acquired during seperate cycles
        img = imread([datadir prefix '/' imgloop(dpc_step).name],...
            'PixelRegion',{[roi(3) roi(4)],[roi(1) roi(2)]});
        img = FilterPixel( img );
        img = img - dark;
        ims(:,:,st,dpc_step) = double( img );
        
        ref = imread([datadir prefix '/' refloop(dpc_step).name],...
            'PixelRegion',{[roi(3) roi(4)],[roi(1) roi(2)]});
        ref = FilterPixel( ref );
        ref = ref - dark;
        
        ims2(:,:,ni*st - 2,dpc_step) = double( img );
        ims2(:,:,ni*st - 1,dpc_step ) = double( ref );
        ims2(:,:,ni*st - 0,dpc_step ) = double( img ) - double( ref );
        
    end
end
fprintf([num2str(num_sett*num_img) ' images loaded, '])
toc
nimplay( ims2(:,:,:,1), 1, [1 2 3], 'step 1' )
nimplay( ims2(:,:,:,3), 1, [1 2 3], 'step 3' )
nimplay( ims2(:,:,:,5), 1, [1 2 3], 'step 5' )

% get settle time
settime = zeros(1,length(prefix_list));
for st = 1:length(prefix_list)
    prefix = prefix_list{st};
    settime(st) = str2double(prefix(20:23));
end

% calculate visibility
vis_map = zeros(roi(4)-roi(3)+1,roi(2)-roi(1)+1,length(prefix_list));
for st = 1:length(prefix_list)
    imsteps = squeeze(ims(:,:,st,:));
    ft_imsteps = fft(imsteps,[],3);
    vis_map(:,:,st) = 2*abs(ft_imsteps(:,:,2))./abs(ft_imsteps(:,:,1));
end

vis_settime = squeeze(mean(mean(vis_map(vis_roi(3):vis_roi(4),...
    vis_roi(1):vis_roi(2),:),1),2));


figure, plot(settime(1:12),vis_settime(1:12),'.-')
hold on
plot(settime(end),vis_settime(end),'*')
xlabel('settle time [ms]')
ylabel('visibility')
ylim([0.06,0.08])

figure, plot(squeeze(ims(vis_roi(3)+5,vis_roi(1)+10,1,:)),'.-')
hold on
plot(squeeze(ims(vis_roi(3)+5,vis_roi(1)+10,2,:)),'.-')
plot(squeeze(ims(vis_roi(3)+5,vis_roi(1)+10,3,:)),'.-')
plot(squeeze(ims(vis_roi(3)+5,vis_roi(1)+10,4,:)),'.-')
plot(squeeze(ims(vis_roi(3)+5,vis_roi(1)+10,5,:)),'.-')
xlabel('phase step')
ylabel('pixel value')
title('stepping curves')
legend({num2str(settime(1:5)')},'Location','eastoutside')

figure, plot(squeeze(ims(vis_roi(3)+5,vis_roi(1):vis_roi(2),1,1)),'.-')
hold on
plot(squeeze(ims(vis_roi(3)+5,vis_roi(1):vis_roi(2),2,1)),'.-')
plot(squeeze(ims(vis_roi(3)+5,vis_roi(1):vis_roi(2),3,1)),'.-')
plot(squeeze(ims(vis_roi(3)+5,vis_roi(1):vis_roi(2),4,1)),'.-')
plot(squeeze(ims(vis_roi(3)+5,vis_roi(1):vis_roi(2),5,1)),'.-')
xlabel('x pixels')
ylabel('pixel value')
title('line profile')
legend({num2str(settime(1:5)')},'Location','eastoutside')

im = ims2(:,:,9,5);
filename = '/asap3/petra3/gpfs/p07/2019/data/11007902/processed/images/ximea_vis04_settime_diff_dpcPos5_settime20ms.png';
im = normat( im );
imwrite( im, filename )

%% Exposure times
datadir = '/asap3/petra3/gpfs/p07/2019/data/11007902/raw/';
prefix_list = { ...
    'ximea_vis07_exptime0001ms_zero', ...
    'ximea_vis07_exptime0002ms_zero', ...
    'ximea_vis07_exptime0003ms_zero', ...
    'ximea_vis07_exptime0005ms_zero', ...
    'ximea_vis07_exptime0010ms_zero', ...
    'ximea_vis07_exptime0015ms_zero', ...
    'ximea_vis07_exptime0020ms_zero', ...
    'ximea_vis07_exptime0025ms_zero', ...
    'ximea_vis07_exptime0030ms_zero', ...
    'ximea_vis08_exptime0001ms_a', ...
    'ximea_vis08_exptime0001ms_b', ...
    'ximea_vis08_exptime0001ms_c', ...
    'ximea_vis08_exptime0001ms_d', ...
    'ximea_vis08_exptime0001ms_e'};
roi = [3750,4050,1350+900,1650+900]; % XIMEA % cropping roi [x1,x2,y1,y2]
vis_roi = [125,175,125+25,175+25]; % XIMEA
num_times = length(prefix_list);
num_img = numel( dir([datadir prefix '/*img*']) );
ims = zeros(roi(4)-roi(3)+1,roi(2)-roi(1)+1,num_times,num_img);
ni = 3;
ims2 = zeros(roi(4)-roi(3)+1,roi(2)-roi(1)+1,num_sett,ni*num_img);
tic
for et = 1:num_times
    prefix = prefix_list{et};
    fprintf( '\n %2u %s:', et,  prefix )
    imgloop = dir([datadir prefix '/*img*']);
    refloop = dir([datadir prefix '/*ref*']);
    darkname = dir([datadir prefix '/*dar*']);
    dark = imread([datadir prefix '/' darkname.name],...
            'PixelRegion',{[roi(3) roi(4)],[roi(1) roi(2)]});
    for dpc_step = 1:num_img
        fprintf( ' %u', dpc_step )
        img = imread([datadir prefix '/' imgloop(dpc_step).name],...
            'PixelRegion',{[roi(3) roi(4)],[roi(1) roi(2)]});
        img = FilterPixel( img );
        img = img - dark;
        ims(:,:,et,dpc_step) = double( img );
                
        ref = imread([datadir prefix '/' refloop(dpc_step).name],...
            'PixelRegion',{[roi(3) roi(4)],[roi(1) roi(2)]});
        ref = FilterPixel( ref );
        ref = ref - dark;
        
        ims2(:,:,ni*et - 2,dpc_step) = double( img );
        ims2(:,:,ni*et - 1,dpc_step) = double( ref );
        ims2(:,:,ni*et - 0,dpc_step ) = double( img ) - double( ref );

    end
end
fprintf([num2str(length(prefix_list)*num_img) ' images loaded, '])
toc

%%
nimplay( cat(2, ims2(:,:,:,1), ims2(:,:,:,2), ims2(:,:,:,3), ims2(:,:,:,4), ims2(:,:,:,5)), 1, [1 2 3], 'step 1' )

% get exposure time
exptime = zeros(1,length(prefix_list));
for et = 1:length(prefix_list)
    prefix = prefix_list{et};
    exptime(et) = str2double(prefix(20:23));
end

figure, subplot(131)
imagesc(ims(:,:,1,1))
title(num2str(exptime(1)))
subplot(132)
imagesc(ims(:,:,3,1))
title(num2str(exptime(3)))
subplot(133)
imagesc(ims(:,:,5,1))
title(num2str(exptime(5)))

% calculate visibility
vis_map = zeros(roi(4)-roi(3)+1,roi(2)-roi(1)+1,length(prefix_list));
for et = 1:length(prefix_list)
    imsteps = squeeze(ims(:,:,et,:));
    ft_imsteps = fft(imsteps,[],3);
    vis_map(:,:,et) = 2*abs(ft_imsteps(:,:,2))./abs(ft_imsteps(:,:,1));
end
for et = 1:length(prefix_list)
    vis_map(:,:,et) = medfilt2(squeeze(vis_map(:,:,et)),[5 5]);
end
vis_exptime = squeeze(mean(mean(vis_map(vis_roi(3):vis_roi(4),...
    vis_roi(1):vis_roi(2),:),1),2));
%%
fig = figure('Name', 'visibility vs exposure time', 'Units', 'Normalized', 'Position', [0.1 0.1 0.5 0.8] );
p = plot(exptime(1:9),vis_exptime(1:9),'x-');
%hold on
%plot(exptime(10:end),vis_exptime(10:end),'*')
xlabel('exposure time [ms]', 'FontSize',font_size )
ylabel('visibility', 'FontSize',font_size )
grid on
ylim([0,0.16])
title( fig.Name, 'FontSize',font_size )
set( p ,'LineWidth', 4, 'MarkerSize', 14 );
ax = gca;
ax.FontSize = font_size;
%legend( 'CdWO_4 100 \mum', 'Location', 'SouthEast' )
axis tight
filename = '/asap3/petra3/gpfs/p07/2019/data/11007902/processed/images/vis_vs_exposure_time.png';
saveas( fig, filename );

%%
figure, plot(squeeze(ims(vis_roi(3)+5,vis_roi(1)+10,1,:)),'.-')
hold on
plot(squeeze(ims(vis_roi(3)+5,vis_roi(1)+10,2,:)),'.-')
plot(squeeze(ims(vis_roi(3)+5,vis_roi(1)+10,3,:)),'.-')
plot(squeeze(ims(vis_roi(3)+5,vis_roi(1)+10,4,:)),'.-')
plot(squeeze(ims(vis_roi(3)+5,vis_roi(1)+10,5,:)),'.-')
xlabel('phase step')
ylabel('pixel value')
title('stepping curves')
legend({num2str(exptime(1:5)')},'Location','eastoutside')

figure, plot(squeeze(ims(vis_roi(3)+5,vis_roi(1):vis_roi(2),1,1)),'.-')
hold on
plot(squeeze(ims(vis_roi(3)+5,vis_roi(1):vis_roi(2),2,1)),'.-')
plot(squeeze(ims(vis_roi(3)+5,vis_roi(1):vis_roi(2),3,1)),'.-')
plot(squeeze(ims(vis_roi(3)+5,vis_roi(1):vis_roi(2),4,1)),'.-')
plot(squeeze(ims(vis_roi(3)+5,vis_roi(1):vis_roi(2),5,1)),'.-')
xlabel('x pixels')
ylabel('pixel value')
title('line profile')
legend({num2str(exptime(1:5)')},'Location','eastoutside')

%% checking aperture size
prefix_list = {'ximea_aperture_08_22ms','ximea_aperture_07_32ms',...
    'ximea_aperture_06_42ms','ximea_aperture_05_55ms',...
    'ximea_aperture_04_82ms','ximea_aperture_03_120ms',...
    'ximea_aperture_02_210ms','ximea_aperture_01_420ms'};

datadir = '/asap3/petra3/gpfs/p07/2019/data/11007902/raw/';

roi = [3450,3750,1350,1650]; % XIMEA % cropping roi [x1,x2,y1,y2]
vis_roi = [100,150,100,150]; % XIMEA

ims = zeros(roi(4)-roi(3)+1,roi(2)-roi(1)+1,length(prefix_list),num_img);
tic
parfor ap = 1:length(prefix_list)
    prefix = prefix_list{ap};
    imgloop = dir([datadir prefix '/*img*']);
    darkname = dir([datadir prefix '/*dar*']);
    dark = imread([datadir prefix '/' darkname.name],...
            'PixelRegion',{[roi(3) roi(4)],[roi(1) roi(2)]});
    for dpc_step = 1:num_img
        this_im = imread([datadir prefix '/' imgloop(dpc_step).name],...
            'PixelRegion',{[roi(3) roi(4)],[roi(1) roi(2)]});
        ims(:,:,ap,dpc_step) = double(this_im)-double(dark);
    end
end
fprintf([num2str(length(prefix_list)*num_img) ' images loaded, '])
toc

% get exposure time
apsize = zeros(1,length(prefix_list));
for ap = 1:length(prefix_list)
    prefix = prefix_list{ap};
    apsize(ap) = str2double(prefix(16:17));
end

exptime = zeros(1,length(prefix_list));
for et = 1:length(prefix_list)
    prefix = prefix_list{et};
    if isnan(str2double(prefix(21)))
        exptime(et) = str2double(prefix(19:20));
    else    
        exptime(et) = str2double(prefix(19:21));
    end
end

meancounts = zeros(1,length(prefix_list));
for ap = 1:length(prefix_list)
    these_ims = squeeze(ims(:,:,ap,:));
    tmp = these_ims(vis_roi(3):vis_roi(4),vis_roi(1):vis_roi(2),:);
    meancounts(ap) = mean(tmp(:));
end

countspers = meancounts./exptime;



vis_range = [1400,2500];
figure, subplot(131)
imagesc(ims(vis_roi(3):vis_roi(4),vis_roi(1):vis_roi(2),1,1),vis_range)
title(num2str(apsize(1)*0.1))
subplot(132)
imagesc(ims(vis_roi(3):vis_roi(4),vis_roi(1):vis_roi(2),round(end/2),1),vis_range)
title(num2str(apsize(round(end/2))*0.1))
subplot(133)
imagesc(ims(vis_roi(3):vis_roi(4),vis_roi(1):vis_roi(2),end,1),vis_range)
title(num2str(apsize(end)*0.1))

% calculate visibility
vis_map = zeros(roi(4)-roi(3)+1,roi(2)-roi(1)+1,length(prefix_list));
for ap = 1:length(prefix_list)
    imsteps = squeeze(ims(:,:,ap,:));
    ft_imsteps = fft(imsteps,[],3);
    vis_map(:,:,ap) = 2*abs(ft_imsteps(:,:,2))./abs(ft_imsteps(:,:,1));
end
for ap = 1:length(prefix_list)
    vis_map(:,:,ap) = medfilt2(squeeze(vis_map(:,:,ap)),[5 5]);
end
vis_apsize = squeeze(mean(mean(vis_map(vis_roi(3):vis_roi(4),...
    vis_roi(1):vis_roi(2),:),1),2));

figure, plot(apsize*0.1,vis_apsize,'.-')
xlabel('aperture size')
ylabel('visibility')
xlim([0,1])
title('cdwo4')

figure, plot(apsize*0.1,vis_apsize.*sqrt(countspers)','.-')
xlabel('aperture')
ylabel('V*sqrt(counts/s)')
title('cdwo4')


figure, yyaxis left
plot(apsize*0.1,vis_apsize,'.-')
ylabel('visibility')
xlabel('aperture size')
yyaxis right
plot(apsize*0.1,countspers,'.-')
ylabel('counts/ms')
xlim([0,1])
title('cdwo4')


%% spare parts below
% make a carpet
talbot_carpet = zeros(roi(2)-roi(1)+1,zpositions);
for zp = 1:zpositions
    talbot_carpet(:,zp) = squeeze(ims(round(end/2),:,zp,1));
end

fun = @(block_struct)  ...
    median(block_struct.data,'all')*ones(size(block_struct.data));
vis_dat = zeros(1,ntomoscans);
for ii = 1:ntomoscans
    vm = vis_map(:,:,ii) ;
    vm2 = blockproc(vm,[5 5],fun);
    vis_dat(ii) = median(vm2,'all');
end


















