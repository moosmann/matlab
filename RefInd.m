% Compute increment of the real part of the refractive index 'delta' from real
% data sets
clear all

water{1} = {'E', 'delta', 'beta'};
water{2} = {'eV', '1', '1'};
water{3} = [20000  5.76620721E-07  3.46201373E-10];
water{4} = [30000  2.56114134E-07  1.06431752E-10];
water{5} = [12000  1.60448508E-06  2.39515319E-09];

fprintf('Water:\nE=%2ukev: %5.3f\nE=%2ukev: %5.3f\nE=%2ukev: %5.3f\n', ...
    water{4}(1)/1000,1e7*water{4}(2),water{3}(1)/1000,1e7*water{3}(2),water{5}(1)/1000,1e7*water{5}(2))
%% Data set. Xenopus 4-cell stage, 40keV, ID19
xen4cell40kev.k = EnergyConverter(40)/2/pi;
xen4cell40kev.pixelsize = 0.75e-6;
xen4cell40kev.hist.tie5p0.full = [-36 31];
xen4cell40kev.hist.tie5p0.roi = [-32 2];
xen4cell40kev.hist.ctf2p5bf0p01.full = [-41 37];
xen4cell40kev.hist.ctf2p5bf0p01.roi = [-38 7];
xen4cell40kev.delta = xen4cell40kev.k*0.745/10000/xen4cell40kev.pixelsize*xen4cell40kev.hist.ctf2p5bf0p01.roi;
[dmin dmax dmean] = Domain(1e7*xen4cell40kev.delta,1,'',0);
fprintf('E=40keV: 10^7*  var(delta)=[%5.2f,%5.2f], max(delta)-min(delta)=%5.2f, mean(delta)=%5.2f\n',dmin,dmax,dmax-dmin,dmean)

%% Data set. Xenopus stage 11, 30keV, BM05
% path: /mnt/tomoraid-LSDF/rci/MI1057-ESRF-BM05-Sept2011/pc/tomo/vol/opt_xeno_stage11_tomo_30keV_20x_500mm_0p30sec_/
xen11E30.k = EnergyConverter(30)/2/pi;
xen11E30.pixelsize = 1.5e-6;
% !! alphaTIE ~= alphaCTF, Reco software was not yet fully consistend
xen11E30.hist.tie5p0.full = [-24 46]; 
xen11E30.hist.tie5p0.roi = [-19 3];
xen11E30.hist.ctf2p5bf0p01.full = [-23 36];
xen11E30.hist.ctf2p5bf0p01.roi = [-19 8];
xen11E30.delta = xen11E30.k*1.5/10000/xen11E30.pixelsize*xen11E30.hist.ctf2p5bf0p01.roi;
[dmin dmax dmean] = Domain(1e7*xen11E30.delta,1,'',0);
fprintf('E=30keV: 10^7*  var(delta)=[%5.2f,%5.2f], max(delta)-min(delta)=%5.2f, mean(delta)=%5.2f\n',dmin,dmax,dmax-dmin,dmean)

%% Data set: Xenopus 4-cell stage, 20keV, ID19
% dat.name = 'Xenopus 4-cell stage';
% dat.path = '/mnt/tomoraid-LSDF/tomo/ESRF_MI1079_ID19_July2011_inlineTomo/Xenopus_4cell/vol/Xenopus_4cell_20keV';
% dat.lambda = EnergyConverter(20);
% dat.pixelsize = 0.75e-6;
% [dat.tie.bins dat.tie.counts] = ReadHist([dat.path '/Histogram of TIElo_alpha5p32.xls']);
% [dat.tie.roi.bins dat.tie.roi.counts] = ReadHist([dat.path '/Histogram of TIElo_alpha5p32_roi.xls']);
% dat.tie.delta = dat.lambda/(2*pi)*0.745/10000/dat.pixelsize*dat.tie.bins;
% dat.tie.roi.delta = dat.lambda/(2*pi)*0.745/10000/dat.pixelsize*dat.tie.roi.bins;
% figure('Name',['Histogram of' dat.name])
% subplot(1,2,1), plot(dat.tie.delta,dat.tie.counts), title('TIE')
% subplot(1,2,2), plot(dat.tie.roi.delta,dat.tie.counts), title('ROI of TIE')
xen4cell20kev.k = EnergyConverter(20)/2/pi;
xen4cell20kev.pixelsize = 0.75e-6;
xen4cell20kev.hist.tie5p32.full = [-65 58];
xen4cell20kev.hist.tie5p32.roi = [-64 -13];
xen4cell20kev.hist.ctf2p5bf0p1.full = [-74 60];
xen4cell20kev.hist.ctf2p5bf0p1.roi = [-72 -10];
xen4cell20kev.delta = xen4cell20kev.k*0.745/10000/xen4cell20kev.pixelsize*xen4cell20kev.hist.ctf2p5bf0p1.roi;
[dmin dmax dmean] = Domain(1e7*xen4cell20kev.delta,1,'',0);
fprintf('E=20keV: 10^7*  var(delta)=[%5.2f,%5.2f], max(delta)-min(delta)=%5.2f, mean(delta)=%5.2f\n',dmin,dmax,dmax-dmin,dmean)

%% Data set. Xenopus stage 10p5, 20keV, ID19
% path: /mnt/tomoraid-LSDF/tomo/ESRF_May2011_Xenopus/vol/Xenopus_10p5Stage_Agar/
xen10p5.k = EnergyConverter(20)/2/pi;
xen10p5.pixelsize = 1.4e-6;
xen10p5.hist.tie5p0.full = [-62 26];
xen10p5.hist.tie5p0.roi = [-51 -10];
xen10p5.hist.ctf2p5bf0p1.full = [-47 17];
xen10p5.hist.ctf2p5bf0p1.roi = [-39 4];
xen10p5.delta = xen10p5.k*1.4/10000/xen10p5.pixelsize*xen10p5.hist.ctf2p5bf0p1.roi;
fprintf('PyHST IMAGE_PIXEL_SIZE_1 value to get correct delta values: %10f\n',10000/(xen10p5.k*1.4/xen10p5.pixelsize))
[dmin dmax dmean] = Domain(1e7*xen10p5.delta,1,'',0);
fprintf('E=20keV: 10^7*  var(delta)=[%5.2f,%5.2f], max(delta)-min(delta)=%5.2f, mean(delta)=%5.2f\n',dmin,dmax,dmax-dmin,dmean)

%% Data set. Xenopus stage 19, 12keV, TopoTomo
% path: /mnt/tomoraid-LSDF/tomo/TopoTomo_100927-SingleDistancePhaseRetrieval/tomo/embryo_st19_50sExT_30cmSD_650proj_360deg_9FFC_9DFI/vol/
xen19.k = EnergyConverter(12)/2/pi;
xen19.pixelsize = 1e-6;
% !! alphaTIE ~= alphaCTF, Reco software was not yet fully consistend
xen19.hist.tie5p0.full = [-147 -196]; %strong tomographic artifacts at corners of image
xen19.hist.tie5p0.roi = [-84 -27];
xen19.hist.ctf2p5bf0p01.roi = [-60 -12];
xen19.delta = xen19.k*1.0/10000/xen19.pixelsize*xen19.hist.ctf2p5bf0p01.roi;
[dmin dmax dmean] = Domain(1e7*xen19.delta,1,'',0);
fprintf('E=12keV: 10^7*  var(delta)=[%5.2f,%5.2f], max(delta)-min(delta)=%5.2f, mean(delta)=%5.2f\n',dmin,dmax,dmax-dmin,dmean)



