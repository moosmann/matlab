function itool(im,varargin)
% Show input image 'im' using imshow with dynamic range [] and
% 'InitialMagnification' 'fit'. Optional input arguments are the name of
% the figure to show (ischar), the slice (isscalar) to show if 'im' is a 3D
% stack, or the ROI (isvector) in the sense im(roi,roi) or im(roi,roi,slice)
%
%Written by Julian Moosmann, last version 2013-10-24

%% Default arguments %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NameOfFigure = inputname(1);
slice = 1;
roi = 0;
dynRange = [];
% Parse input arguments
for nn = 1:numel(varargin)  
    arg = varargin{nn};    
    if ischar(arg)        
        NameOfFigure = arg;
       % fprintf('Name of figure: %s\n',NameOfFigure);
    end
    
    if isnumeric(arg)
        switch numel(arg)
            case 1
                slice = arg;
                % fprintf('Show slice of input stack: %u\n',slice);
            case 2
                dynRange = arg;
            otherwise
                roi = arg;
                % fprintf('Show ROI: %u:%u\n',roi(1),roi(end));
        end               
    end
end
imroislice = @(im) im(roi,roi,slice);
imslice = @(im) im(:,:,slice);

%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
im = gather( im );
if roi(1) > 0
    h = imtool( imroislice( squeeze( im ) ),dynRange,'InitialMagnification','fit');
    set(h,'units','normalized','outerposition',[0.05 0.05 0.95 0.95]);
    set(h,'Name',sprintf('Image Tool: %s.  ROI: %u:%u.  Input size: %ux%ux%u %',NameOfFigure,roi(1),roi(end),size(im),slice))
else
    h = imtool( imslice( squeeze(im) ),dynRange,'InitialMagnification','fit');
    set(h,'units','normalized','outerposition',[0.05 0.05 0.95 0.95]);
    set(h,'Name',sprintf('Image Tool: %s.  Input size: %ux%ux%u %',NameOfFigure,size(im),slice))
end
