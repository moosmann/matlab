function im = SetDynRange(im,DynRange)
% Set dynamic range of input image 'im' to be within [DynRange(1)
% DynRange(2)].
%
% Written by Julian Moosmann, 2013-12-09
%
% im = SetDynRange(im,DynRange) 

%% Default Arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    DynRange = [0 1];
end

%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

im(im < DynRange(1)) = DynRange(1);
im(im > DynRange(2)) = DynRange(2);

