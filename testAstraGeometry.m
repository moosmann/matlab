% Test volume and projection geometries

vol_geom = astra_create_vol_geom(256, 256);
proj_geom = astra_create_proj_geom('parallel', 1.0, 256, linspace2(0,pi,180));

% As before, create a sinogram from a phantom
P = phantom(256);
[sinogram_id, sinogram] = astra_create_sino_gpu(P, proj_geom, vol_geom);
%figure(1); imshow(P, []);
figure(1); imshow(sinogram, []);

proj_geom = astra_create_proj_geom('parallel', 1/2.0, 512, linspace2(0,pi,180));
[sinogram_id, sinogram] = astra_create_sino_gpu(P, proj_geom, vol_geom);
figure(2); imshow(sinogram, []);
