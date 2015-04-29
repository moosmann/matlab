% Load original image. 
%load  woman
% Generate noisy image.
%x = X + 15*randn(size(X));

filepath = '/mnt/LSDF/tomo/APS_32ID-C_LifeCellImaging_GUP31523_2012-10-13/savedLocally/phase/Oct11_02-48_wildtype_stage11p5_34p5keV_0700mm_20ms_0834proj_scantime50s_deadtime10p5min_20ms_open_40ms_close/3DstackProc_FiltSino_FDcor_tie_regPar2p50_noMeanSub/tomo00/phase_0001.edf';
x = pmedfread(filepath)';
[dimx dimy] = size(x);
xcrop = ceil(dimx/8):floor(7/8*dimx);
ycrop = ceil(dimy/8):floor(7/8*dimy);
x = x(xcrop,ycrop);

filepath = '/mnt/LSDF/tomo/APS_32ID-C_LifeCellImaging_GUP31523_2012-10-13/savedLocally/vol/Oct11_02-48_wildtype_stage11p5_34p5keV_0700mm_20ms_0834proj_scantime50s_deadtime10p5min_20ms_open_40ms_close/3DstackProc_FiltSino_FDcor_tie_regPar2p50_noMeanSub/';
filepath = [filepath 'tomo00_slice0651.tif'];
x = imread(filepath);
% Find default values. In this case fixed form threshold
% is used with estimation of level noise, thresholding
% mode is soft and the approximation coefficients are 
% kept.
[thr,sorh,keepapp] = ddencmp('den','wv',x);

% thr is equal to estimated_sigma*sqrt(log(prod(size(X))))

% De-noise image using global thresholding option.
xd = wdencmp('gbl',x,'sym4',2,thr,sorh,keepapp);
%% Plots.
itool(x)
itool(xd)
itool(x-xd)
% colormap(pink(255)), sm = size(map,1);
% subplot(221), image(wcodemat(X,sm)), title('Original Image')
% subplot(222), image(wcodemat(x,sm)), title('Noisy Image')
% subplot(223), image(wcodemat(xd,sm)), title('denoised Image')