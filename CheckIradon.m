clear all
ca
InputPath = '/mnt/tomoraid-LSDF/tomo/APS_32ID-C_LifeCellImaging_GUP31523_2012-10-13/savedLocally/phase/Oct11_00-30_wildtype_stage11p0_34p5keV_0700mm_15ms_0834proj_scantime50s_deadtime08min_20ms_open_40ms_close/3DstackProc_FiltSino_FDcor_tie_regPar2p50_noMeanSub/tomo00_tif/';
OutputPath = '/mnt/tomoraid-LSDF/users/moosmann/sino/';
adjMeth = 2;
%% Make sinogram
tic
% sino: HorizontalDimOfProjection x NumberOfProjections (x
% slicePerpendicularToRotationAxis)
sino = MakeSino(601,InputPath,'phase_*')'; 
tsino = toc;
%% Adjust dimension of sinogram
tic;
[dimHor NumProj] = size(sino);
sinoCen = dimHor/2;
RotAxisPos = 725.50;
pixshift = abs(round(2*(sinoCen-RotAxisPos)));
switch adjMeth
    case {1,'crop'}
        % Cropping
        if round(RotAxisPos) < round(sinoCen)
            sino     = sino(1:end-pixshift,:);
        elseif round(RotAxisPos) > round(sinoCen)
            sino     = sino(1+pixshift:end,:);
        end
    case {2,'pad'}
        % Padding
        dimHorPad = 2048;
        if round(RotAxisPos) < round(sinoCen)
            sino = padarray(sino,[ceil((dimHorPad-dimHor)/2-pixshift/2) 0],'symmetric','pre');
            sino = padarray(sino,[floor((dimHorPad-dimHor)/2+pixshift/2) 0],'symmetric','post');
            
        elseif round(RotAxisPos) > round(sinoCen)
            sino = padarray(sino,[ceil( (dimHorPad-dimHor)/2 +pixshift/2) 0],'symmetric','pre');
            sino = padarray(sino,[floor((dimHorPad-dimHor)/2 -pixshift/2) 0],'symmetric','post');
        end
end
dimHorNew = size(sino,1);
tadj = toc;
%% Inverse radon trafo
theta = 180/(NumProj-1);
interp = 'linear';
filter = 'Ram-Lak';
freq_scaling = 1;
%output_size = 2*floor(dimHorNew/(2*sqrt(2)));% MATLAB'S default
output_size = 1*dimHorNew;
fprintf('DimHor: %g, DimHorNew: %g, output_size: %g\n',dimHor,dimHorNew,output_size);
tic;
vol = iradon(sino,theta,interp,filter,freq_scaling,output_size);
ttomo = toc;
%% 
fprintf('Time for making sinogram: %g s, adjusting sinogram dimension: %g s, inverse radon transform: %g s.\n',tsino,ttomo);
itool(vol)