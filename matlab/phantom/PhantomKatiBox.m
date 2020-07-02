dataPath = '/home/jmoosmann/data/phantom/Kati_box/';

s = load([dataPath 'Kati_box.mat']);

k = s.Kati_box;
khu = s.Kati_box_HU;

s = load([dataPath 'Sinogram_greyscale_nsdv1procent.mat']);
s = s.sinogram;

fprintf('\nphantom\n size: %u, %u, %u',size(k))
fprintf('\n     min: %g, max: %g', min(k(:)), max(k(:)))
fprintf('\n HU: min: %g, max: %g', min(khu(:)), max(khu(:)))
fprintf('\nsinogram\n size: %u, %u, %u',size(s))
fprintf('\n min: %g, max: %g', min(s(:)), max(s(:)))
fprintf('\n')


