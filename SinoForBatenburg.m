if 0
    OutputPath = '/mnt/tomoraid-LSDF/users/moosmann/Batenburg/LowQualityDataSet/';
    
    InputPath = '/mnt/tomoraid-LSDF/tomo/APS_2BM_LifeCellImaging_GUP28266/lifecellimaging/phase/wildtype_30keV_10min_deadtime_25tomo_stage11p0_upwards_620mm_015ms/int_filtLineSectionMedWidthH063V001_noCropping/tie_alpha2p50_padNo_filtRing/tomo01/';
    
    MakeSino([441:460 671:692],InputPath,'phase_*',[OutputPath 'sino'],'tif');
end
if 1
    OutputPath = '/mnt/tomoraid-LSDF/users/moosmann/Batenburg/MediumQualityDataSet/';
    InputPath = '/mnt/tomoraid-LSDF/tomo/APS_32ID-C_LifeCellImaging_GUP31523_2012-10-13/savedLocally/phase/Oct11_00-30_wildtype_stage11p0_34p5keV_0700mm_15ms_0834proj_scantime50s_deadtime08min_20ms_open_40ms_close/3DstackProc_FiltSino_FDcor_tie_regPar2p50_noMeanSub/tomo00/';
    VolPath = '/mnt/tomoraid-LSDF/tomo/APS_32ID-C_LifeCellImaging_GUP31523_2012-10-13/savedLocally/vol/Oct11_00-30_wildtype_stage11p0_34p5keV_0700mm_15ms_0834proj_scantime50s_deadtime08min_20ms_open_40ms_close/3DstackProc_FiltSino_FDcor_tie_regPar2p50_noMeanSub';
    MakeSino(601:620,InputPath,'phase_*',[OutputPath 'sino'],'tif');
end