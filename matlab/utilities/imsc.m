function imsc(varargin)
% imagesc using gray colormap.
%
% Written by Julian Moosmann, 2015-03-16. Modified: 2016-12-01

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imagesc(varargin{:})
colormap(gray)
