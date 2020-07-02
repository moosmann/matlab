% Size in microns of cell from tomographic reconstrucion of stage 12
eggsize = (1900-120)*0.75; % microns
% During gastrulation cell moving from the outer edge to the interior. This
% movement is about a quarter of the egg.
movlen = 1/4*eggsize; % microns
movtim = 3*3600; % seconds
% Print.
format short
fprintf('\nsize of egg from tomogram: %.3f millimeter\n',eggsize/1000)
fprintf('cell movement during gastrulation: %.3f microns\n',movlen)
fprintf('cell movement velocity: %.6f microns/s\n',movlen/movtim)