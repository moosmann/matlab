%% dpc scripts
% data is in '/asap3/petra3/gpfs/p07/2019/data/11007902/raw/'

datadir = '/asap3/petra3/gpfs/p07/2019/data/11007902/raw/test_vis_01/';

a = dir([datadir '*img*']);
nims = numel( a );
sy = 3840;
sx = 5120;
ims = zeros(sy,sx,nims);
tic
parfor ii = 1:nims
    tmp = imread([datadir filesep a(ii).name ]);
    ims(:,:,ii) = tmp;
end
toc

stepping_curve = squeeze(ims(2319,2215,:));
figure, plot(stepping_curve,'.-')


roi = [1750,2050,2050,2350]; % [x1,x2,y1,y2]

figure, imagesc(ims(roi(3):roi(4),roi(1):roi(2),1))

figure, plot(squeeze(ims(roi(3),2050:2150,1)))



crop_im = squeeze(ims(roi(3):roi(4),roi(1):roi(2),1));
ft_crop_im = fftshift(fft2(crop_im));

figure, imagesc(log(abs(ft_crop_im))), axis equal


crp_im = ims(roi(3):roi(4),roi(1):roi(2),1:5);
ft_crp_im = fft(crp_im,[],3);

vis_map = abs(ft_crp_im(:,:,2))./abs(ft_crp_im(:,:,1));

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
nps = numel(tmp2); % number of phase steps

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

ims = zeros(roi(4)-roi(3)+1,roi(2)-roi(1)+1,zpositions,nps);
tic
parfor zp = 1:zpositions
    innerloop = dir([datadir outerloop(zp).name '/*img*']);
    darkname = dir([datadir outerloop(zp).name '/*dar*']);
    dark = imread([datadir outerloop(zp).name '/' darkname.name],...
            'PixelRegion',{[roi(3) roi(4)],[roi(1) roi(2)]});
    for ps = 1:nps
        tmp = imread([datadir outerloop(zp).name '/' innerloop(ps).name],...
            'PixelRegion',{[roi(3) roi(4)],[roi(1) roi(2)]});
        ims(:,:,zp,ps) = double(tmp-dark);
    end
end
fprintf([num2str(zpositions*nps) ' images loaded, '])
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

z = 800;%:10:1000;
tiltangles = 32:1:44;

zpositions = length(z);
tilts = length(tiltangles);
ntomoscans = numel(outerloop);
nps = numel(dir([datadir outerloop(1).name '/*img*'])); % number of phase steps

% load data
refname = dir([datadir outerloop(1).name '/*ref*']);
tmp = imread([datadir outerloop(1).name '/' refname(1).name]);
tmp_crop = tmp(roi(3):roi(4),roi(1):roi(2));
[sy,sx,sz] = size(tmp);
roi = [3750,4050,1000,1300]; % XIMEA % cropping roi [x1,x2,y1,y2]
vis_roi = [125,175,125+25,175+25]; % XIMEA
figure, imagesc(tmp), axis equal
hold on, rectangle('Position',[roi(1),roi(3),roi(2)-roi(1),roi(2)-roi(1)],'EdgeColor','r')
title('crop roi')
figure, imagesc(tmp_crop), axis equal
hold on, rectangle('Position',[vis_roi(1),vis_roi(3),vis_roi(2)-vis_roi(1),...
    vis_roi(2)-vis_roi(1)],'EdgeColor','r')
title('vis roi')

ims = zeros(vis_roi(4)-vis_roi(3)+1,vis_roi(2)-vis_roi(1)+1,ntomoscans);
tic
for ii = 1:ntomoscans
    innerloop = dir([datadir outerloop(ii).name '/*img*']);
    darkname = dir([datadir outerloop(ii).name '/*dar*']);
    dark = imread([datadir outerloop(ii).name '/' darkname.name],...
            'PixelRegion',{[roi(3) roi(4)],[roi(1) roi(2)]});
    for ps = 1:nps
        tmp = imread([datadir outerloop(ii).name '/' innerloop(ps).name],...
            'PixelRegion',{[roi(3) roi(4)],[roi(1) roi(2)]});
        tmp = double(tmp-dark);
        ims(:,:,ii,ps) = tmp(vis_roi(3):vis_roi(4),vis_roi(1):vis_roi(2));
    end
end
fprintf([num2str(ntomoscans*nps) ' images loaded, '])
toc

vis_map = zeros(vis_roi(4)-vis_roi(3)+1,vis_roi(2)-vis_roi(1)+1,ntomoscans);
for ii = 1:ntomoscans
    these_ims = squeeze(ims(:,:,ii,:));
    ft_ims = fft(these_ims,[],3);
    vm = 2*abs(ft_ims(:,:,2))./abs(ft_ims(:,:,1));
    vis_map(:,:,ii) = medfilt2(vm,[3 3]);
end

vis_dat = squeeze(mean(mean(vis_map,1),2));
figure, plot(tiltangles,vis_dat,'.-')
xlabel('tilt angle [deg]')
ylabel('visibility')

vis_dat = reshape(vis_dat,[zpositions,tilts]);



figure, imagesc(z,tiltangles,vis_dat)
xlabel('z [mm]')
ylabel('tilt angle [deg]')
title(prefix,'Interpreter','none')


%% checking settle times
prefix_list = {'ximea_vis04_settime0000ms','ximea_vis04_settime0010ms',...
    'ximea_vis04_settime0020ms','ximea_vis04_settime0030ms','ximea_vis04_settime0040ms',...
    'ximea_vis04_settime0050ms','ximea_vis04_settime0070ms',...
    'ximea_vis04_settime0080ms','ximea_vis04_settime0090ms','ximea_vis04_settime0100ms',...
    'ximea_vis04_settime0150ms','ximea_vis04_settime0200ms',...
    'ximea_vis06_settime0020ms_zero'};

datadir = '/asap3/petra3/gpfs/p07/2019/data/11007902/raw/';

roi = [3750,4050,1350+900,1650+900]; % XIMEA % cropping roi [x1,x2,y1,y2]
vis_roi = [125,175,125+25,175+25]; % XIMEA

ims = zeros(roi(4)-roi(3)+1,roi(2)-roi(1)+1,length(prefix_list),nps);
tic
parfor st = 1:length(prefix_list)
    prefix = prefix_list{st};
    innerloop = dir([datadir prefix '/*img*']);
    darkname = dir([datadir prefix '/*dar*']);
    dark = imread([datadir prefix '/' darkname.name],...
            'PixelRegion',{[roi(3) roi(4)],[roi(1) roi(2)]});
    for ps = 1:nps
        tmp = imread([datadir prefix '/' innerloop(ps).name],...
            'PixelRegion',{[roi(3) roi(4)],[roi(1) roi(2)]});
        ims(:,:,st,ps) = double(tmp-dark);
    end
end
fprintf([num2str(length(prefix_list)*nps) ' images loaded, '])
toc

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

%% checking exposure times
prefix_list = {'ximea_vis07_exptime0001ms_zero','ximea_vis07_exptime0002ms_zero'...
    'ximea_vis07_exptime0003ms_zero','ximea_vis07_exptime0005ms_zero',...
    'ximea_vis07_exptime0010ms_zero','ximea_vis07_exptime0015ms_zero',...
    'ximea_vis07_exptime0020ms_zero','ximea_vis07_exptime0025ms_zero',...
    'ximea_vis07_exptime0030ms_zero',...
    'ximea_vis08_exptime0001ms_a','ximea_vis08_exptime0001ms_b',...
    'ximea_vis08_exptime0001ms_c','ximea_vis08_exptime0001ms_d',...
    'ximea_vis08_exptime0001ms_e'};

datadir = '/asap3/petra3/gpfs/p07/2019/data/11007902/raw/';

roi = [3750,4050,1350+900,1650+900]; % XIMEA % cropping roi [x1,x2,y1,y2]
vis_roi = [125,175,125+25,175+25]; % XIMEA

ims = zeros(roi(4)-roi(3)+1,roi(2)-roi(1)+1,length(prefix_list),nps);
tic
parfor et = 1:length(prefix_list)
    prefix = prefix_list{et};
    innerloop = dir([datadir prefix '/*img*']);
    darkname = dir([datadir prefix '/*dar*']);
    dark = imread([datadir prefix '/' darkname.name],...
            'PixelRegion',{[roi(3) roi(4)],[roi(1) roi(2)]});
    for ps = 1:nps
        this_im = imread([datadir prefix '/' innerloop(ps).name],...
            'PixelRegion',{[roi(3) roi(4)],[roi(1) roi(2)]});
        ims(:,:,et,ps) = double(this_im)-double(dark);
    end
end
fprintf([num2str(length(prefix_list)*nps) ' images loaded, '])
toc

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


figure, plot(exptime(1:9),vis_exptime(1:9),'.-')
hold on
plot(exptime(10:end),vis_exptime(10:end),'*')
xlabel('exposure time [ms]')
ylabel('visibility')


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

ims = zeros(roi(4)-roi(3)+1,roi(2)-roi(1)+1,length(prefix_list),nps);
tic
parfor ap = 1:length(prefix_list)
    prefix = prefix_list{ap};
    innerloop = dir([datadir prefix '/*img*']);
    darkname = dir([datadir prefix '/*dar*']);
    dark = imread([datadir prefix '/' darkname.name],...
            'PixelRegion',{[roi(3) roi(4)],[roi(1) roi(2)]});
    for ps = 1:nps
        this_im = imread([datadir prefix '/' innerloop(ps).name],...
            'PixelRegion',{[roi(3) roi(4)],[roi(1) roi(2)]});
        ims(:,:,ap,ps) = double(this_im)-double(dark);
    end
end
fprintf([num2str(length(prefix_list)*nps) ' images loaded, '])
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

figure, plot(apsize*0.1,vis_apsize.*sqrt(countspers)','.-')
xlabel('aperture')
ylabel('V*sqrt(counts/s)')

figure, yyaxis left
plot(apsize*0.1,vis_apsize,'.-')
ylabel('visibility')
xlabel('aperture size')
yyaxis right
plot(apsize*0.1,countspers,'.-')
ylabel('counts/ms')
xlim([0,1])


%% spare parts below
% make a carpet
talbot_carpet = zeros(roi(2)-roi(1)+1,zpositions);
for zp = 1:zpositions
    talbot_carpet(:,zp) = squeeze(ims(round(end/2),:,zp,1));
end



















