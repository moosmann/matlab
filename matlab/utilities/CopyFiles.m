function CopyFiles(DoCopy)

if nargin < 1
    DoCopy = 0;
end
%% Parameters
% Path = '/mnt/tomoraid-LSDF/tomo/APS_2BM_LifeCellImaging_GUP28266/lifecellimaging/data/wildtype_30keV_10min_deadtime_20tomo_stage14p0_upwards_620mm_025ms';
% SourceFileRange = 18*1420+1300+60+(1:60);
% TargetFileRange = 19*1420+1300+60+(1:60);
% Path = '/mnt/tomoraid-LSDF/tomo/APS_2BM_LifeCellImaging_GUP28266/lifecellimaging/data/wildtype_30keV_10min_deadtime_10tomo_stage12p5_upwards_620mm_050ms';
% SourceFileRange = 8*1420+1300+60+(1:60);
% TargetFileRange = 9*1420+1300+60+(1:60);
% Path = '/mnt/tomoraid-LSDF/tomo/APS_2BM_LifeCellImaging_GUP28266/lifecellimaging/data/wildtype_30keV_05min_deadtime_08tomo_stage11p0_upwards_620mm_050ms';
% TomoNum = 8;
% SourceFileRange = (TomoNum-2)*1420+1300+60+(1:60);
% TargetFileRange = (TomoNum-1)*1420+1300+60+(1:60);
Path = '/mnt/tomoraid-LSDF/tomo/APS_2BM_LifeCellImaging_GUP28266/lifecellimaging/data/wildtype_30keV_05min_deadtime_08tomo_stage11p0_upwards_620mm_050ms';
TomoNum = 8;
SourceFileRange = (TomoNum-2)*1420+1300+60+(1:60);
TargetFileRange = (TomoNum-1)*1420+1300+60+(1:60);
%% Body
fprintf('Source file numbers: %u - %u, %u files\n',SourceFileRange(1),SourceFileRange(end),numel(SourceFileRange))
fprintf('Target file numbers: %u - %u, %u files\n',TargetFileRange(1),TargetFileRange(end),numel(TargetFileRange))
for nn = 1:numel(TargetFileRange)
    SourceName = sprintf('%s/proj_%05u.tif',Path,SourceFileRange(nn));
    TargetName = sprintf('%s/proj_%05u.tif',Path,TargetFileRange(nn));
    fprintf('Copy process number: %2u. Exists source file: %u. Exists target file: %u\n%s\n%s\n',nn,exist(SourceName,'file'),exist(TargetName,'file'),SourceName,TargetName)
    if DoCopy
        copyfile(SourceName,TargetName,'f');
    end
end