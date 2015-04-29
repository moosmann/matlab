% %% Script is written to remove salt and pepper noise from the data
% %% Adviced for Fotron camera
% %% Written by V.Altapova, August 2010, first version

function imfiltered=filter_im1(img,threshold)
   
    if nargin<2, threshold = 1 ; end;
    
    img_median = medfilt2(img,[3 3]);
    
    R=img./img_median;
    [m n]=size(img);
    imfiltered=img;
    s=0;
    for i=1:m
        for j=1:n
            if img(i,j)>=threshold
               s=s+1; 
                imfiltered(i,j)=img_median(i,j);
                
            end
        end
    end
    %disp('Number of corrected pixels'); s
    
