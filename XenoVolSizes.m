clear all

set(1).name = ('wildtype_30keV_10min_deadtime_25tomo_stage11p0_upwards_620mm_015ms');
set(1).numvol = 9;
set(1).volsizeorg = 6.4;
set(1).volsizereduced = 1.1;
set(1).volsizeunit = 'GB';
fprintf('total size of all used volumes: %g GB\n',set(1).numvol*set(1).volsizeorg)
fprintf('reduced size of all used volumes: %g GB\n',set(1).numvol*set(1).volsizereduced)