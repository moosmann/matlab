function [sino, expStr,sinoStr] = getSino(Experiment,SubsetIndex,Slice)

if nargin < 3
    Slice = 0;
end

switch lower(Experiment)
    case {'xeno4','xenopus4','xenopus4cell','xeno4cell'}
        expStr = 'xenopus4cell';
        % Subsets of different filterings and phase retrievals
        sinoName{1} = 'sino_int_phase_tie_regPar2p50_tif';
        sinoName{2} = 'sino_int_phase_quasi_regPar2p50_binFilt0p100_tif';
        sinoName{3} = 'sino_int_phase_quasiNew_regPar2p50_binFilt0p100_tif';
        sinoName{4} = 'sino_int_phase_filtSino_tie_regPar2p50_tif';
        sinoName{5} = 'sino_int_phase_filtSino_quasi_regPar2p50_binFilt0p100_tif';
        sinoName{6} = 'sino_int_phase_filtSinoTrans_tie_regPar2p50_tif';        
        sinoName{7} = 'sino_int_phase_filtSinoTrans_quasi_regPar2p50_binFilt0p100_tif';
        sinoName{8} = 'sino_int_phase_filtSinoTrans_quasiNew_regPar2p50_binFilt0p100_tif';
        sinoName{9} = 'sino_int_phase_filtSinoTrans_ctfHalfSine_regPar2p50_tif';
        sinoName{10} = 'sino_int_phase_filtSinoBoth_tie_regPar2p50_tif';
        sinoName{11} = 'sino_int_phase_filtSinoBoth_quasi_regPar2p50_binFilt0p100_tif';
        sinoName{12} = 'sino_int_tif';
        sinoName{13} = 'sino_int_phase_quasi_regPar2p00_binFilt0p100_tif';
        
        % Choose correct parentPath
        if regexp(sinoName{SubsetIndex},'phase')
            parentPath = '/export/scratch1/moosmann/ESRF_MI1079_ID19_July2011_inlineTomo/phase/Xenopus_4cell_20keV/';
        else 
            parentPath = '/export/scratch1/moosmann/ESRF_MI1079_ID19_July2011_inlineTomo/int/Xenopus_4cell_20keV/';                        
        end
        sinoStr = sinoName{SubsetIndex};
        sinoPath = sprintf('%s%s/',parentPath,sinoStr);
        filenameCell = FilenameCell(sinoPath);
        if Slice
            filename = sprintf('%s%s',sinoPath,filenameCell{Slice});
        else
            filename = sprintf('%s%s',sinoPath,filenameCell{round(numel(filenameCell)/2)});
        end
        sino = imread(filename);
        % Last four projections are for control
        sino = sino(1:1:end-4,:);
        % Crop sinogram such that rot axis positon is in center
        rotpos = 1067;
        sino = RotAxisSymmetricCropping(sino,rotpos);
        % Delete write format from folder name
        sinoStr(end-3:end) = [];
end

