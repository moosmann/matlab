phaStr  = '/mnt/tomoraid3/tomo/APS_2BM_LifeCellImaging_GUP28266/phase/wildtype_05min_deadtime_05tomo_stage16p0_upwards_620mm_050ms_30p0keV/TIElo_alpha4p02/tomo01/phase_0001.edf';
intStr = '/mnt/tomoraid3/tomo/APS_2BM_LifeCellImaging_GUP28266/int/wildtype_05min_deadtime_05tomo_stage16p0_upwards_620mm_050ms_30p0keV/tomo01/int_0001.edf';
pha  = pmedfread(phaStr)';
flatStr = '/mnt/tomoraid3/tomo/APS_2BM_LifeCellImaging_GUP28266/data/wildtype_05min_deadtime_05tomo_stage16p0_upwards_620mm_050ms_30p0keV/proj_01400.tif';
flat = double(imread(flatStr));
flat = FilterHotPixel(flat,0.03);
%ishow(pha),ishow(flat)
%I = imread('rice.png');
I = flat;
figure(1), imshow(I,[],'InitialMagnification','fit')
background = imopen(I,strel('disk',20));
figure(2),imshow(background,[],'InitialMagnification','fit')
%figure, surf(double(background(1:8:end,1:8:end))),zlim([0 255]);
set(gca,'ydir','reverse');
I2 = I - background;
figure(3),imshow(I2,[],'InitialMagnification','fit')
I3 = imadjust(I2);
%ishow(I3);