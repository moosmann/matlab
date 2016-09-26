% /mnt/tomoraid-LSDF/tomo/APS_2BM_LifeCellImaging_GUP28266/lifecellimaging/vol/wildtype_30keV_10min_deadtime_06tomo_stage11p0_upwards_620mm_050ms

numProj = 1200;
overallAngle = 180;
angleInc =  0.15;
% Number of projection where projection of volume 6 corresponds to
% projection #1 of volume 5
numIm = 52 - 1;
shiftVol6toVol5 = -angleInc*numIm;
disp(shiftVol6toVol5)


% /mnt/tomoraid-LSDF/tomo/APS_2BM_LifeCellImaging_GUP28266/lifecellimaging/vol
numIm = 162 - 1;
shiftVol6toVol5 = -angleInc*numIm;
disp(shiftVol6toVol5)

% /mnt/tomoraid-LSDF/tomo/APS_2BM_LifeCellImaging_GUP28266/lifecellimaging/vol/
% wildtype_30keV_10min_deadtime_25tomo_stage11p0_upwards_620mm_015ms/int_filtLineSectionMedWidthH063V001_noCropping/tomo10

numIm = 321 - 1;
shiftVol6toVol5 = -angleInc*numIm;
disp(shiftVol6toVol5)

 % Suppl Video 5b (4b), stage 12.5, folder: wildtype_30keV_10min_deadtime_10tomo_stage12p5_upwards_620mm_050ms 
 % offset of 1st projections wrt 2nd projections 
 numIm = 54 - 1;
shiftVol6toVol5 = angleInc*numIm;
disp(shiftVol6toVol5)