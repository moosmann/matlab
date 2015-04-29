function [vol,angles]= astra_vol(sino,filter_type,angularIncRad,projIndexInc)


tic
[NumProj,NumPix,NumSlices] = size(sino);

if nargin < 2
    filter_type = 'ram-lak';
end
if nargin <3
    angularIncRad = 2*pi/1599;
end
if nargin < 4
    projIndexInc = round(1599/NumProj);
end

angles = projIndexInc * angularIncRad*((1:NumProj)-1);

for nn = NumSlices:-1:1
    %PrintNum(nn,30)
    % Creat volume and projection geometry
    vol_geom = astra_create_vol_geom([NumPix]);
    %proj_geom = astra_create_proj_geom('parallel',1.0,NumPix,2*pi/1599*linspace(0,2*pi,NumProj));
    proj_geom = astra_create_proj_geom('parallel',1.0,NumPix,angles);
    % Create the sinogram data object
    
    sino_id = astra_mex_data2d('create', '-sino', proj_geom, SubtractMean(squeeze(sino(:,:,nn))));
    % Create reconstruction
    [vol_id, vol(:,:,nn)] = test_astra_fbp_cuda(proj_geom,vol_geom,sino_id,filter_type);    
    astra_mex_data2d('delete', vol_id);
    astra_mex_data2d('delete', sino_id);
end

fprintf('\n Elapsed time: %g s\n',toc);

%pha = ifftn( PhaseFilter3D('tie',size(vol),[20 0.945 0.75e-6],2.5,0.1,'single' ) .* fftn(vol))
%dx = round(1.1*(size(vol,1)-size(vol,1)/sqrt(2))/2);x=dx:size(vol,1)-dx;nimplay(pha(x,x,:))