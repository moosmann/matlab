function im = FilterMedian(im,Neigborhood_Ver_Hor,Padding)
% Perform median filtering of input matrix.
%
% im: matrix.
% Neigborhood_Ver_Hor: 1x2-vector. default: [3 3]
% Padding: string or scalar. default: 'symmetric'. other: 'zeros',
% 'indexed'
%
% Written by Julian Moosmann, last version: 2013-11-07
%
% im = FilterMedian(im,Neigborhood_Ver_Hor,Padding)

%% Defaults
if nargin < 2;
    Neigborhood_Ver_Hor = [3 3];
end
if nargin < 3
    Padding = 'symmetric';
end

%% Main
im = medfilt2(im,Neigborhood_Ver_Hor,Padding);