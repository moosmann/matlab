function imsc(varargin)
% imagesc using gray colormap and bilinear interpolation.

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imagesc(varargin{:},'Interpolation','bilinear')
colormap(gray)
