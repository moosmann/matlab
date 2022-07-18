function ishow(im,varargin)
% Show input image 'im' using imshow with dynamic range [] and
% 'InitialMagnification' 'fit'. Optional input arguments are the name of
% the figure to show (ischar), the slice (isnumeric & numel==1) to show if 'im' is a 3D
% stack, or the ROI (isnumeric & numel==2) in the sense im(roi,roi) or im(roi,roi,slice)
%
%Written by Julian Moosmann, last version 2019-01-30

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default arguments
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
im = squeeze(im);
rotate_im = size(im,1) > size(im,2) ;
rot = @(im) rot90( im, rotate_im );
if rotate_im
    fprintf( '\n Image is rotated to fit to the screen!' )
end
%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if roi(1) > 0
    im = rot( imroislice( squeeze( im ) ) );
    figure('Name',sprintf('%s.  ROI: %u:%u.  Input size: %ux%ux%u %',NameOfFigure,roi(1),roi(end),size(im),slice))
    %set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    imshow( im, dynRange,'InitialMagnification','fit', 'XData',[1 size(im,2)], 'YData',[1 size(im,1)] );
    colorbar;
    axis equal tight
else
    im = rot( imslice( squeeze( im ) ) );
    figure('Name',sprintf('%s.  Input size: %ux%ux%u %',NameOfFigure,size(im),slice))
    %set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    imshow( im, dynRange,'InitialMagnification','fit','XData',[1 size(im,2)], 'YData',[1 size(im,1)] );
    colorbar;    
    axis equal tight
end
