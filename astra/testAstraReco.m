clear all
%% Read sinogram
%sino=Readstack('/export/scratch1/moosmann/art/Batenburg/MediumQualityDataSet/sino',100);
%rotpos = 725.5;
parentPath = '/export/scratch1/moosmann/ESRF_MI1079_ID19_July2011_inlineTomo/phase/Xenopus_4cell_20keV/';
sinoPhase{1} = 'sino_int_phase_filtSino_tie_regPar2p50_tif';
sinoPhase{2} = 'sino_int_phase_quasiNew_regPar2p50_binFilt0p100_tif';
sinoPhase{3} = 'sino_int_phase_quasi_regPar2p50_binFilt0p100_tif';
sinoPhase{4} = 'sino_int_phase_filtSino_quasi_regPar2p50_binFilt0p100_tif';
sinoPhase{5} = 'sino_int_phase_filtSinoTrans_quasi_regPar2p50_binFilt0p100_tif';
sinoPhase{6} = 'sino_int_phase_filtSinoTrans_quasiNew_regPar2p50_binFilt0p100_tif';
sinoPhase{7} = 'sino_int_phase_filtSinoTrans_ctfHalfSine_regPar2p50_tif';
sinoPhase{8} = 'sino_int_phase_filtSinoBoth_tie_regPar2p50_tif';
sinoPhase{9} = 'sino_int_phase_filtSinoBoth_quasi_regPar2p50_binFilt0p100_tif';

sinoPhaseIndex = [5];

for nn = numel(sinoPhaseIndex):-1:1
    sinoPath = [parentPath sinoPhase{sinoPhaseIndex(nn)} '/'];
    filenameCell = FilenameCell(sinoPath);
    filename = sprintf('%s%s',sinoPath,filenameCell{round(numel(filenameCell)/2)});
    rotpos = 1067;
    sino = imread(filename);
    sino = sino(1:1:end-4,:);
    % Crop sinogram such that rot axis positon is in center
    if rotpos < size(sino,2)/2
        sino = sino(:,1:round(2*rotpos));
    else
        horWidth = size(sino,2)-rotpos;
        sino = sino(:,end-round(2*horWidth):end);
    end
    sino = sino(1:4:end,1:2:end);
    dimHor = size(sino,2);
    sinos(:,:,nn) = sino;
    %sino = RotAxisSymmetricPadding(sino,rotpos,0);
    
    % Filter phase map
    sino = RemoveLowFreq(sino(:,256:end-256),9,[1 18]);
    
    % Retrieved phase maps have zero mean, but that doesn't imply zero mean of
    % the detector lines of the sinogram
    
    % Creat volume and projection geometry
    %vol_geom = astra_create_vol_geom([size(sino,2),size(sino,2)]);
    %dimHorVol = 2*round(dimHor/sqrt(2)/2);
    dimHorVol = dimHor;
    vol_geom = astra_create_vol_geom([dimHorVol,dimHorVol]);
    proj_geom = astra_create_proj_geom('parallel',1.0,size(sino,2),linspace2(0,2*pi,size(sino,1)));
    
    
    % We now re-create the sinogram data object as we would do when loading
    % an external sinogram
    sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, sino);
    
    % Create a data object for the reconstruction
    rec_id = astra_mex_data2d('create', '-vol', vol_geom);
    
    % Set up the parameters for a reconstruction algorithm using the GPU
    %cfg = astra_struct('FBP_CUDA');
    cfg = astra_struct('SIRT_CUDA');
    cfg.ReconstructionDataId = rec_id;
    cfg.ProjectionDataId = sinogram_id;
    % Use GPU #1 for the reconstruction. (The default is #0.)
    cfg.Options.GPUindex = 1;
    
    % Available algorithms:
    % SIRT_CUDA, SART_CUDA, EM_CUDA, FBP_CUDA (see the FBP sample)
    
    
    % Create the algorithm object from the configuration structure
    alg_id = astra_mex_algorithm('create', cfg);
    
    % Run 150 iterations of the algorithm
    tic;
    astra_mex_algorithm('iterate', alg_id, 100);
    toc;
    
    % Get the result
    rec(:,:,nn) = astra_mex_data2d('get', rec_id);
    rec_par_string = sprintf('%s_proj%04u' ,cfg.type,size(sino,1));
    if dimHorVol < dimHor
        x = 1:dimHorVol;
    else
        x=round(dimHorVol*(1-1/sqrt(2))):round(dimHorVol/sqrt(2));
    end
    %itool(rec(x,x),1,rec_par_string);
    %itool(rec,1,rec_par_string);
end
if size(rec,3) > 1
    nimplay(rec(x,x,:),0)
else
    itool(rec(x,x),1,rec_par_string)
end
% Clean up. Note that GPU memory is tied up in the algorithm object,
% and main RAM in the data objects.
astra_mex_algorithm('delete', alg_id);
astra_mex_data2d('delete', rec_id);
astra_mex_data2d('delete', sinogram_id);
