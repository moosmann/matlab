% %% Script is written to remove salt and pepper noise from the data
% %% Adviced for Fotron camera
% %% Written by V.Altapova, August 2010, modified by J.Moosmann, Jan 2011
% Final version, so far. Compare to filter_im.m

function [imfiltered,s,img_median,m]=filter_im3(img,threshold,img_median)

    if nargin<2, threshold = 1 ; end
    if nargin<3, img_median = medfilt2(img,[3 3],'symmetric');end

    R=img./img_median;
    imfiltered=img;
    imfiltered(R>threshold)=img_median(R>threshold);
    s=numel(imfiltered(R>threshold));
    m=zeros(size(img));m(R>threshold)=1;

    %disp('Number of corrected pixels'); s

