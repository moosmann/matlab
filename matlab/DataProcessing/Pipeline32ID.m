function Pipeline32ID(useGPU,doDarkFieldCorrection,LineFiltering,doFilterSino,doMeanSubstraction,doNormalizeSlices)
% Script which (interactively) loops over the scans of the 'Comparative
% Life-cell imaging' experiment at 32ID@APS, Chigaco, in October 2012

%% Default arguments
if nargin < 1
    useGPU = 0;
end
if nargin < 2
    doDarkFieldCorrection = 1;
end
if nargin < 3
    LineFiltering = [0 0 0];
end
if nargin < 4
    doFilterSino = 1;
end
if nargin < 5
    doMeanSubstraction = 0;
end
if nargin < 6
    doNormalizeSlices = 0;
end
showFigures = 0;

%% Logging of output
[~, hostname] = unix('echo $HOSTNAME');
eval(sprintf('diary %sMatlabOutputLog_%s_%s.txt','/mnt/tomoraid-LSDF/users/moosmann/matlabLog/',hostname,datestr(now,'yyyy-mm-dd_HH-MM-SS')));
fprintf('\nSCRITP: ''Pipeline32ID''\nSTART: %s\n',datestr(now,'yyyy.mm.dd  HH:MM:SS'));
fprintf('HOSTNAME: %s',hostname);
%% Parameters
ParentPath = '/mnt/tomoraid-LSDF/tomo/APS_32ID-C_LifeCellImaging_GUP31523_2012-10-13/savedLocally/';
out = ChooseDatSet([ParentPath 'data']);
TomoSetsToProcess   = input(sprintf('TOMOGRAMS TO RECONSTRUCT [0..], DEFAULT=ALL: '));
if useGPU
    matlabpool close force local
end
if ~doDarkFieldCorrection
    fprintf('\n\nNO DARK FIELD CORRECTION WILL BE DONE\n\n')
end
%% Loop over different stages/embryos
for nn = 1:numel(out)
    DataSet = out(nn).name;
    %% Start loop over time-lapsed tomograms
    %DataProc32ID(TomoSetsToProcess,useGPU,doDarkFieldCorrection,LineFiltering,showFigures,'DataSet',DataSet,'ParentPath',ParentPath,'PipeGPUtag',1);
    DataProc32IDstack(TomoSetsToProcess,useGPU,doDarkFieldCorrection,LineFiltering,doFilterSino,doMeanSubstraction,doNormalizeSlices,showFigures,'DataSet',DataSet,'ParentPath',ParentPath,'PipeGPUtag',useGPU)
end
%% Stop matlabpool and logging
if useGPU
    matlabpool close force local
end
diary off