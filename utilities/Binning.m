function im = Binning(im, bin)
% 2 x 2 or 4 x 4 binning of 2-D image. Images is cropped before binning
% such that mod(size(im), bin) = 0. 
%
% im: 2D image to bin
% bin: bin size, 2 or 4. Default: 2
%
% Written by Julian Moosmann.
% Last modification 2016-12-06
%
% im = Binning(im, bin)

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    bin = 2;
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if bin == 0 || bin == 1    
    return    

elseif bin == 2
    % crop to even number of pixels    
    im = im(1:size(im,1)-mod(size(im,1),2),1:size(im,2)-mod(size(im,2),2)); 

    % bin
    im = im(1:2:end,1:2:end) + im(2:2:end,1:2:end) + im(1:2:end,2:2:end) + im(2:2:end,2:2:end);
    
elseif bin == 4
    % crop to even number of pixels    
    im = im(1:size(im,1)-mod(size(im,1),4),1:size(im,2)-mod(size(im,2),4)); 

    % bin
    im = im(1:2:end,1:2:end) + im(2:2:end,1:2:end) + im(1:2:end,2:2:end) + im(2:2:end,2:2:end);
    im = im(1:2:end,1:2:end) + im(2:2:end,1:2:end) + im(1:2:end,2:2:end) + im(2:2:end,2:2:end);        

else
    error('Bin size %g not implemented', bin)
end
