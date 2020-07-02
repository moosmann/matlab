%InputPath = '/mnt/LSDF/tomo/APS_32ID-C_LifeCellImaging_GUP31523_2012-10-13/savedLocally/phase/Oct11_00-30_wildtype_stage11p0_34p5keV_0700mm_15ms_0834proj_scantime50s_deadtime08min_20ms_open_40ms_close/3DstackProc_FiltSino_FDcor_tie_regPar2p50_noMeanSub/tomo00';
InputPath = '/mnt/tomoraid-LSDF/tomo/APS_2BM_LifeCellImaging_GUP28266/lifecellimaging/phase/wildtype_30keV_10min_deadtime_25tomo_stage11p0_upwards_620mm_015ms/int_filtLineSectionMedWidthH063V001_noCropping_filtSino/TIElo_alpha4p02/tomo01';
OutputPath = [InputPath '_tif'];

CheckAndMakePath(OutputPath);
FileStruct = dir([InputPath '/phase*']);
parfor nn=1:numel(FileStruct)
    filepath = [InputPath '/' FileStruct(nn).name];
    im = pmedf_read_jm(filepath)';
    filepath = [OutputPath '/' FileStruct(nn).name];
    filepath = [filepath(1:end-3) 'tif'];
    write32bitTIF(filepath,im);
end