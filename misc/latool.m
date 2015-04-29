function latool(im,varargin)
% Take log of abs of input image 'im', and call itool.
%
% Written by Julian Moosmann, 2013-10-14, last version 2013-10-24
%
% latool(im,varargin)

itool(log(1 + abs(im)),varargin)