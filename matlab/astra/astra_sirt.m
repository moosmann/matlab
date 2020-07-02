% Directory
parfolder = '/mntdirect/_data_visitor/ls2395/id19/';
name = '20th_01h03m58_Stage21_Gap22p5_sequ1sample_2_';'19th_04h06m12_Stage24_Gap22p5_sequ1sample_1_';
folder = [parfolder name];
% Projection paramtere
numProj = 500;
angles = ((1:500) - 1) / 500 * pi;
% ROI
%x = 1:1600;
x = 500:1100;
y = 1:1280;
numIter = 150;
dimPad = 270;

% Flat field
if 0
    ref = double(Readstack(folder, 1,'ref*_0000.edf'));
    flat = median(ref, 3);
    fprintf('\n')
    fprintf('\n flat: %u x %u', size(flat))
    flatfilt = FilterPixel(flat);
    flatroi = flatfilt(x, y);
    fprintf('\n flatroi: %u x %u', size(flatroi))
end

% Sino
if 1
    pf = PhaseFilter('qphalfsine', size(flatroi), [22.16 5.43 1.6e-6], 2.5, 0.1, 'double')';
    for nn = numProj:-1:1
        filename = sprintf('%s/%s%04u.edf', folder, name, nn);
        im = ( FilterPixel( double(edfread(filename, y, x))' ) ./ flatroi )';
        sino3d(nn, :, :) = im;
        %sino3d(nn, :, :) = ifft2( pf .* fft2( im ) );
        if 0
            im = edfread(filename, y, x)'./flat(x, y);
            imagesc(im)
            colormap gray
            pause(0.01)
        end
    end
    fprintf('\n sino3d: %u x %u x %u', size(sino3d))
    sino = SubtractMean(squeeze(sino3d(:,:,300)));
end

% Loop plot over projections
if 0
    for nn = 1:numProj
        imagesc(squeeze(sino3d(nn, :, :))')
        colormap gray
        axis equal tight
        pause(0.1)
    end
end


recfbp = astra_make_fbp(sino, angles);

% Configuration of algorithmus
dimHor = size(sino,2);
volDim = dimHor;
dimHorPad = dimHor + 2 * dimPad;

% Create geometries
proj_geom = astra_create_proj_geom('parallel', 1.0, dimHor, angles);
proj_geom_pad = astra_create_proj_geom('parallel', 1.0, dimHorPad, angles);
vol_geom  = astra_create_vol_geom(volDim, volDim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% row sum of system matrix via FP of unit volume

% store volume of ones
vol_id = astra_mex_data2d('create', '-vol', vol_geom, 1);

% Create forward projection
sino_id = astra_mex_data2d('create', '-sino', proj_geom, 0);
cfg = astra_struct('FP_CUDA');
cfg.ProjectionDataId = sino_id;
cfg.VolumeDataId = vol_id;
cfg.option.GPUindex = 0;
fp_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('run', fp_id);
rowsum = 1./astra_mex_data2d('get', sino_id);
aclear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% column sum of system matrix via BP of unit sinogram

% Store sino of ones
sino_id = astra_mex_data2d('create', '-sino', proj_geom, 1);

% Create backprojection
reco_id = astra_mex_data2d('create', '-vol', vol_geom, 0);
cfg = astra_struct('BP_CUDA');
cfg.ProjectionDataId = sino_id;
cfg.ReconstructionDataId = reco_id;
cfg.option.GPUindex = 0;
bp_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('run', bp_id);
colsum = 1./astra_mex_data2d('get', reco_id);
aclear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIRT
figure('Name', 'SIRT iteration')
for nn = 1:numIter
    % Initialize volume
    if nn == 1
        rec = zeros(volDim, volDim);
    end
    
    % 1: Forward projection of current reconstruction
    % Store current volume
    vol_id = astra_mex_data2d('create', '-vol', vol_geom, rec);
    % Create forward projection
    sino_id = astra_mex_data2d('create', '-sino', proj_geom, 0);
    cfg = astra_struct('FP_CUDA');
    cfg.ProjectionDataId = sino_id;
    cfg.VolumeDataId = vol_id;
    cfg.option.GPUindex = 0;
    fp_id = astra_mex_algorithm('create', cfg);
    astra_mex_algorithm('run', fp_id);
    p = astra_mex_data2d('get', sino_id);
    aclear
    
    % 2: Subtract result from data
    p = sino - p;
    
    % 3: Multiply with rowsum of system matrix
    p = rowsum .* p;
    
    % 4: Back projection of result
    % Store result
    sino_id = astra_mex_data2d('create', '-sino', proj_geom_pad, padarray(p, [0 dimPad], 'replicate', 'both'));
    % Create backprojection
    reco_id = astra_mex_data2d('create', '-vol', vol_geom, 0);
    cfg = astra_struct('BP_CUDA');
    cfg.ProjectionDataId = sino_id;
    cfg.ReconstructionDataId = reco_id;
    cfg.option.GPUindex = 0;
    bp_id = astra_mex_algorithm('create', cfg);
    astra_mex_algorithm('run', bp_id);
    
    % 5: Multiply with colsum of system matrix and update
    rec = rec + colsum .* astra_mex_data2d('get', reco_id);
    aclear
    
    subplot(1,1,1)
    imagesc(rec)
    title(sprintf('iteration: %u', nn))
    colormap gray
    axis equal tight
    pause(0.1)
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Garbage disposal
astra_mex_data2d('clear');
astra_mex_algorithm('clear');



% Plot
if 0
    figure(1)
    subplot(1,2,1)
    imagesc(recfbp)
    colormap gray
    axis equal tight
    subplot(1,2,2)
    imagesc(recsirt)
    colormap gray
    axis equal tight
end