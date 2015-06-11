PhasePath = '/mnt/tomoraid-LSDF/tomo/APS_2BM_LifeCellImaging_GUP28266/lifecellimaging/phase/wildtype_30keV_10min_deadtime_25tomo_stage11p0_upwards_620mm_015ms/TIElo_alpha4p02/tomo01';
TarPath = '/mnt/tomoraid-LSDF/tomo/APS_2BM_LifeCellImaging_GUP28266/lifecellimaging/test';

for nn=1:1200
    source = sprintf([PhasePath '/phase_%04u.edf'],nn);
    target = sprintf([TarPath '/phase_%04u.edf'],1201-nn);
    copyfile(source,target)
    %fprintf('%s\n%s\n',source,target)
end