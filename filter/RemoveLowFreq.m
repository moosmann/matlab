function im = RemoveLowFreq(im,sigma,hsize,BoundaryOption)
% Deprecated. Name changed to 'FilterLowFreq'. See 'FilterLowFreq' for
% details.

%% Defaults arguments
if nargin < 2
   sigma = 40;
end
if nargin < 3
   hsize = sigma*[3 3];
end
if nargin<4
    BoundaryOption = 'replicate';
end
%% MAIN
im = FilterLowFreq(im,sigma,hsize,BoundaryOption);