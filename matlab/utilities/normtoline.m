function im = normtoline(im)
% Normalize 2D image line-wise.
%
% writtten by Julian Moosmann, 2016-09-26
% last version: 2016-09-26


%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
im = squeeze(im);

im = im ./ repmat( mean( im, 1), [size(im, 1) 1] );
