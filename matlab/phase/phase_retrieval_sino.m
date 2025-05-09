clear all
ca
%% Parameter
path_scan = '/asap3/petra3/gpfs/p07/2024/data/11020415/processed/itaw008_cet_337a_hi1405_dl';
path_sino = [path_scan '/trans02_180'];
path_sino_phase = [path_scan '/trans02_180_phase'];
%path_reco_phase = [path_scan '/reco_phase'];
take_neg_log = 1;
edp = [37e3 1 1.27e-6];
regpar = 1;
bf = 0;
fc = 1;
prec = 'single';
padding = 2;
visual_output = 1;
window_state = 'minimized';'normal';'maximized';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
t = toc;
%CheckAndMakePath(path_reco_phase)
CheckAndMakePath(path_sino_phase)
fprintf('\nPhase retrieval for sinograms')
fprintf('\n sinogram path: %s',path_sino)
d = dir([path_sino,'/*tif']);
num_sino = numel(d); % height
n = round(num_sino/2);
fn = [d(n).folder filesep d(n).name];
imfo = imfinfo(fn);
sino = imread(fn,'Info',imfo);
num_pixel = imfo.Width; % width
num_angle = imfo.Height; % projections
folder = d(1).folder;
name = {d(:).name};
proj_shape = [num_sino,num_pixel];
dt = class(sino);
fprintf('\n sinogram size: [angles pixels] = [%u %u]',num_angle,num_pixel)
% Plot sinogram
if visual_output
    if exist('h0' , 'var' ) && isvalid( h0 )
        figure(h0)
    else
        h0 = figure('Name', 'sinogram: transmission and phase', 'WindowState', window_state );
    end
    subplot(1,2,1)
    imsc(sino);
    title(sprintf('transmission sinogram'))
    colorbar
    axis equal tight
    xticks('auto'),yticks('auto')
    drawnow
end

fprintf('\n subregion:')
% horizontal pixel row subregion
p0 = 1;
dp = 1;
p1 = num_pixel;
p = p0:dp:p1;
num_p = numel(p);
fprintf('\n  horizontal pixels: %u of %u',num_p,num_pixel)
% vertical pixel column subregion
y0 = 1;
dy = 1;
y1 = num_sino;
y = y0:dy:y1;
num_y = numel(y);
fprintf('\n  horizontal pixels (sinograms): %u of %u',num_y,num_sino)
% phase filter
[pf,pfs] = PhaseFilter('tie',padding * [num_y,num_p],edp,regpar,bf,fc,prec);

%% Angle slab loop
angle_slab_width = 50; % angle/proj slab width
num_angle_slabs = ceil(num_angle/angle_slab_width);
for ai = 1:num_angle_slabs
    ts = toc;
    % angular subregion
    a0 = 1 + (ai-1)*(angle_slab_width);
    da = 1;
    a1 = (a0-1) + angle_slab_width;
    a1 = min([a1,num_angle]);
    a = a0:da:a1;
    num_a = numel(a);
    fprintf('\n  angles (projections): %u of %u',num_a,num_angle)
    fprintf('\n  angle loop %3u: [start increment end] = [%u %u %u]',ai,a0,da,a1)
    % preallocation
    proj = zeros([num_y,num_p,num_a],dt);
    fprintf('\n projection slab size: [height pixel angles] = [%u %u %u]',size(proj))

    %% Read projection from sinograms
    fprintf('\n Reading projections from sinograms')
    t = toc;
    fprintf('\n  ')
    c = 0;
    for yi = 1:num_y
        fn = [folder filesep name{y(yi)}];
        s = imread(fn,'tif','Info',imfo,'PixelRegion',{[a0 da a1],[p0 dp p1]});
        s = NegLog(s,take_neg_log);
        proj(yi,:,:) = shiftdim(s',-1);
        if mod(yi-1,100) == 0
            fprintf(' %u',yi);
            c = c + 1;
            if mod(c,20) == 0
                fprintf('\n  ')
            end
        end
    end
    fprintf('\n  duration: %.0f s = %.1f min',toc-t, (toc-t)/60)
    if visual_output
        h1 = figure('Name', 'projection and phase', 'WindowState', window_state );
        subplot(1,2,1)
        imsc(proj(:,:,1));
        title(sprintf('projection'))
        colorbar
        axis equal tight
        xticks('auto'),yticks('auto')
        drawnow
    end

    %% Phase retrieval
    t = toc;
    fprintf('\n Phase retrieval')
    pha = padarray( proj, (padding-1) * [num_y,num_p], 'symmetric', 'post' );
    pha = fft2(pha);
    pha = pf .* pha ;
    pha = ifft2(pha);
    pha = -real(pha);
    pha = pha(1:num_y,1:num_p,:);
    fprintf('\n  duration: %.0f s = %.1f min',toc-t, (toc-t)/60)
    if visual_output
        if exist('h1' , 'var' ) && isvalid( h1 )
            figure(h1)
        else
            h1 = figure('Name', 'projection and phase', 'WindowState', window_state );
        end
        subplot(1,2,2)
        imsc(pha(:,:,1));
        title(sprintf('phase'))
        c = gray;
        c = flipud(c);
        colormap(c)
        colorbar
        axis equal tight
        xticks('auto'),yticks('auto')
        drawnow
    end

    %% Save sinograms
    t = toc;
    fprintf('\n Write phase sinograms')
    start = [a0 1];
    count = [num_a,num_p];
    stride = [1 1];
    c = 0;
    fprintf('\n Writing h5 phase sinograms from phase maps')
    fprintf('\n ')
    for h = 1:num_y
        fnh = sprintf('%s/sino_phase_%06u.h5',path_sino_phase,h);
        if isequal(ai,1)
            if exist(fnh,'file')
                s = h5info(fnh,'/sino_phase');
                if ~isequal(s.Dataspace.Size,[num_angle,num_p])
                    delete(fnh);
                    h5create(fnh, '/sino_phase', [num_angle,num_p], 'Datatype', dt);
                end
            else
                h5create(fnh, '/sino_phase', [num_angle,num_p], 'Datatype', dt);
            end
        end
        s = shiftdim(squeeze(pha(h,:,:)),1);
        h5write(fnh, '/sino_phase',s,start,count,stride);
        if mod(h-1,100) == 0
            fprintf(' %u',h);
            c = c + 1;
            if mod(c,20) == 0
                fprintf('\n  ')
            end
        end
    end
    fprintf('\n  duration: %.0f s = %.1f min',toc-t, (toc-t)/60)
    fprintf('\n angle slab %u of %u finished in %.0f = %.1f',ai,num_angle_slabs,toc-ts,(toc-ts)/60)
end
% Plot phase sinogram
hd = dir([path_sino_phase '/*.h5']);
hnum = numel(hd);
hn = round(hnum/2);
if visual_output
    if exist('h0' , 'var' ) && isvalid( h0 )
        figure(h0)
    else
        h0 = figure('Name', 'sinogram: transmission and phase', 'WindowState', window_state );
    end
    subplot(1,2,2)
    fnh = [hd(hn).folder filesep hd(hn).name];
    hsino = h5read(fnh,'/sino_phase');
    imsc(hsino);
    title(sprintf('phase sinogram'))
    colorbar
    axis equal tight
    xticks('auto'),yticks('auto')
    drawnow
end
fprintf('\n  total duration: %.0f s = %.1f min',toc-t0, (toc-t0)/60)
fprintf('\n Finished\n')