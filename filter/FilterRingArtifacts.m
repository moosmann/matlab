clear all;ca
%% Read volumes slice
VolDir = '/mnt/tomoraid-LSDF/tomo/APS_2BM_LifeCellImaging_GUP28266/lifecellimaging/test/vol/int_filtLineSectionMedWidthH065V001_defaultCropping/wildtype_30keV_10min_deadtime_25tomo_stage11p0_upwards_620mm_015ms/';
VolName = 'TIElo_alpha4p02_tomo01.vol';
dimx = 1024;
dimy = 1024;
dimz = 1008;
VolDim = [dimx dimy dimz];
DataType = 'float32';
SubVolPos = [1 1 ceil(dimz/2)];
SubVolDim = [dimy dimx 1];
im = double(mexVolRead([VolDir VolName],VolDim,DataType,SubVolPos,SubVolDim));
%% Transformation to polar coordinates
RotAxisPos = 489.499989;
if dimx/2 > RotAxisPos
    xmin = RotAxisPos;
else
    xmin = dimx-RotAxisPos;
end
if dimy/2 > RotAxisPos
    ymin = RotAxisPos;
else
    ymin = dimy-RotAxisPos;
end
rmax = floor(min([xmin ymin])-1);
rsteps = 1*floor(rmax);
radius = 0:rmax/rsteps:rmax*(1-1/rsteps);
thetasteps = 360;%2*pi*rmax;
theta = 0:2*pi/thetasteps:2*pi*(1-1/thetasteps);

tic;
for rr = length(radius):-1:1
    for tt = length(theta):-1:1
        impc(rr,tt) = im(round(RotAxisPos+radius(rr)*sin(theta(tt))),round(RotAxisPos+radius(rr)*cos(theta(tt))));
    end
end
tpolar = toc;
fprintf('Polar coordinate transform of [DimX DimY] = [%u %u] with rotation axis at %g to [DimRadius DimAngle] = [%u %u] in %gs.\n',...
    VolDim(1),VolDim(2),RotAxisPos,size(impc),tpolar)
tic;
%% Radial median filtering
MedFiltRadius = 15;
MedFiltTheta  = 1;
impcmf = medfilt2(impc,[MedFiltRadius MedFiltTheta],'symmetric');
tmedfilt = toc;
fprintf('Median filtering with mask [%u %u] done in %gs.\n',MedFiltRadius,MedFiltTheta,tmedfilt)
% ishow(impc)
% ishow(impcmf)
% ishow(impc-impcmf)
%% Inverse transform to cartesian coordinates
xmax = cos(pi/4)*radius(end);
for xx = -xmax+1:xmax
    for yy = -xmax+1:xmax
        imf(xx,yy) = impc(round(sqrt(xx^2+yy^2)),round(atan(yy/xx)));
    end
end