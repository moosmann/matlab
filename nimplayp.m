function nimplayp(imstack,normGlobal,permuteOrder)
% Call nimplay with predefined permuation order. See function nimplay for
% details.
%
% Written by Julian Moosmann, last version 2013-10-18

if nargin < 2
    normGlobal = 1;
end
if nargin < 3
    permuteOrder = [3 2 1];
end

fprintf('Rearrange array dimensions by permutation vector %s\n',mat2str(permuteOrder));
%% Call nimplay
nimplay(imstack,normGlobal,permuteOrder);