%vector containing the row number that should be read
RowsToRead = 1:1008;
%string of the absolute or relative input path. Default is the current folder
InputDir = '/mnt/tomoraid-LSDF/tomo/APS_2BM_LifeCellImaging_GUP28266/lifecellimaging/int/wildtype_30keV_10min_deadtime_25tomo_stage11p0_upwards_620mm_015ms/int_filtLineSectionMedWidthH063V001_noCropping/tomo01';
%string containg the naming pattern
ImageNamePattern = [];
% String of the output path where the sinograms will be stored
OutputDirSino = [];
OutputDirImages = '/mnt/tomoraid-LSDF/tomo/APS_2BM_LifeCellImaging_GUP28266/lifecellimaging/int/wildtype_30keV_10min_deadtime_25tomo_stage11p0_upwards_620mm_015ms/int_filtLineSectionMedWidthH063V001_noCropping_filtSino/tomo01';
ProjectionsToRead = 1:1200;
HotPixThresh = 0;
%%
ScriptFilterSino(RowsToRead,InputDir,ImageNamePattern,OutputDirSino,OutputDirImages,ProjectionsToRead,HotPixThresh);