function nimplayp(imstack,renorm_slicewise,permuteOrder)
% Call nimplay with predefined permuation order. See function nimplay for
% details.
%
% Written by Julian Moosmann, last version 2013-10-18, 2017-04-27

if nargin < 2
    renorm_slicewise = 1;
end
if nargin < 3
    permuteOrder = [2 1 3];
end

fprintf('Rearrange array dimensions by permutation vector %s\n',mat2str(permuteOrder));
%% Call nimplay
nimplay(imstack,renorm_slicewise,permuteOrder);