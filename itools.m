function itools(varargin)
%Show images given as input with Magnification ('fit') and Contrast ([])
%adjusted.

slice = 1;
%% Body
for nn=1:numel(varargin)
    im = squeeze(varargin{nn});
    h = imtool(im(:,:,slice),[],'InitialMagnification','fit');
    set(h,'Name',['Image Tool: ' inputname(nn)]);
end
