InputPath = '/mnt/tomoraid-LSDF/users/moosmann/CWI_DATA/ESRF_MI1079_ID19_July2011_inlineTomo/vol/Xenopus_4cell_20keV/tomo_intProjInc1_fbp__1962x1962x0100/';

if ~exist('stack','var')
    stack = Readstack(InputPath);
end


%% Pad volume
[dimx dimy dimz] = size(stack);
stack = padarray(stack,[dimx/2 dimy/2 dimz/2],'symmetric','both');

%% FT
stack = fftn(stack);

%% PhaseFilter
stack = stack .* PhaseFilter3D('qp',size(stack),[20 0.945 0.745e-6],3,0.1,'single');

%% iFT
stack = ifftn(stack);

stack = real(stack);

%% Crop
stack = stack(dimx/2 + (1:dimx),dimy/2 +(1:dimy),dimz/2 +(1:dimz));


