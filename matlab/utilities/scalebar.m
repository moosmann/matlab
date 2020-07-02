% -------------------------------------------------------------------------
% scalebar
% -------------------------------------------------------------------------
% Michael Hughes, February 2017
% michael.robert.hughes@gmail.com                     www.mike-hughes.org
% -------------------------------------------------------------------------
% Adds a scale bar to the pixel data of an image. By adding the scale
% bar directly to the image pixel data, rather than as an overlay in the
% figure window, the scale bar can then be saved along with the image in a 
% bitmap file format without altering the image data.
%
% The scalebar can also be labelled with text. This is done via a hack 
% to allow text to be converted to a bitmap. It requires opening and 
% closing a figure window which cannot be hidden without producing unstable
% results. Inevitablly, for small images the text will appear poor
% quality and pixelated.
% -------------------------------------------------------------------------
% Usage:
%
%      imageOut = scalebar(image, barWidth, barLength, label, varargin)
%
%   image       : image to add scalebar to. Either monochrome or colour
%   barWidth    : width of scale bar in pixels (y dimension in image)
%   barLength   : length of scale bar in pixels (x dimension in image) -
%                  this is the distance that provides the 'scale'.
%   label       : text label to place next to the bar. Set as '' for no
%                 label
%   varargin    : any of the following optional name-pair paramaeters:-
%
%       position       : position of scale bar. Allowed: 'bottom-right' 
%                        [Default],'top', 'bottom', 'left', 'right',  
%                        'top-left', 'top-right','bottom-left' or 'exact'. 
%                        If 'exact' use barX and barY to specify position.
%       barX           : left position of scale bar in pixels if 'position'
%                        is 'exact'.
%       barY           : top position of scale bar in pixels if 'position' 
%                        is 'exact'.
%       colour         : colour of scale bar, vector of length 3 specifying 
%                        pixel values in the R, G and B planes. If image is 
%                        mono then only the R values is used. Default is 
%                        [255, 255, 255].
%       textBoxPosition: position of text relative to scale bar. 'top', 
%                        'bottom' (default), 'left' or 'right'.
%       textBoxPadding : space in pixels between text and scale bar. 
%                        Default is 10
%       fontName       : font to use for text label. Default is 'Arial'.
%       fontSize       : size of font to use for text label. Default is 15.
%       fontWeight     : 'bold' or 'normal' (default).
%       fontAngle      : 'italic' or 'normal' (default).
% -------------------------------------------------------------------------
 
function image = scalebar(image, barWidth, barLength, label, varargin)

% Parse optional name-pairs
options = struct('position' , 'bottom-right', 'textBoxPosition', 'bottom', 'colour', [255 255 255], 'edgeMargin', 10, ... 
        'textBoxPadding', 10, 'fontName' , 'Arial', 'fontSize' , 15, 'fontWeight', 'normal', 'fontAngle', 'normal', ...
        'barX', 1, 'barY', 1);
    
optionNames = fieldnames(options);
for pair = reshape(varargin,2,[])
    inpName = pair{1};
    if any(strcmp(inpName, optionNames))
        options.(inpName) = pair{2};
    else
        error('Invalid parameter.');
    end
end

colour = options.colour;
textBoxPosition = options.textBoxPosition;
position = options.position;
edgeMargin = options.edgeMargin;
textBoxPadding = options.textBoxPadding;
fontName = options.fontName;
fontSize = options.fontSize;
fontWeight = options.fontWeight;
fontAngle = options.fontAngle;
barX = options.barX;
barY= options.barY;

% For use of the 'text' function, we need the text colour to be
% in the range 0 to 1. So work out how to scale the text colour and
% keep a record so we can convert it to the pixel value we really want
% at the end when we copy the text onto the image
colourScale = 1./max(colour);
if isinf(colourScale)
    colourScale = 1;
end    
textcolour = colour * colourScale;

% Check if we need to add a text label
isLabel = ~strcmp(label,'');

% Pull out dimensions of image
width = size(image,2);
height = size(image,1);
centreX = size(image,2) / 2;
centreY = size(image,1) / 2;

% If there is a text label, figure out size of text box
if isLabel

   % Draw a blank image (value 255) in a figure
   h = figure;
   blank = zeros(size(image,1),size(image,2)); blank(:) = 255; 
   imshow(blank);
  
   textTestX = 20;
   textTestY = 20;
   
   % Write the text in the centre of the image
   t = text(textTestX, textTestY, label, 'FontName' , fontName, 'FontSize', fontSize, 'FontWeight', fontWeight, 'FontAngle', fontAngle, 'VerticalAlignment', 'top', 'Margin', 1);
 
   % Extract the frame containing the imag
   fr = getframe(gca,[10,10,size(blank,2)-10,size(blank,1)-10]);

   % Get rid of the text
   delete(t);
   
   % Extract the pixel data
   c = fr.cdata(:,:,1);

   % Find all the pixels which contain text
   [textPixelsY,textPixelsX] = find(c~=255);
      
   % Work out the dimension of the rectangle enclosing the text
   left = min(textPixelsX);
   right = max(textPixelsX);
   top = min(textPixelsY);
   bottom = max(textPixelsY);
   textOffsetX = left - textTestX;
   textOffsetY = top - textTestY;
   textBoxHeight = bottom - top;
   textBoxWidth = right-left;

   % Set the margins to be applied to the bar, depending on the size
   % of the text box
   rightTextMargin = max(0,textBoxWidth - barLength) / 2;  
   leftTextMargin = max(0,textBoxWidth - barLength) / 2;
   topTextMargin = max(0,textBoxHeight - barWidth) / 2 ;
   bottomTextMargin = max(0,textBoxHeight - barWidth) /  2;
   
   switch(textBoxPosition)
        case 'top'
            topTextMargin = textBoxHeight + textBoxPadding;
            bottomTextMargin = 0;
        case 'bottom'
            bottomTextMargin = textBoxHeight + textBoxPadding;
            topTextMargin = 0;       
        case 'left'
            leftTextMargin = textBoxWidth + textBoxPadding;
            rightTextMargin = 0;
        case 'right'
            rightTextMargin = textBoxWidth + textBoxPadding;
            leftTextMargin = 0;
    end
    
else % If we have no text box, don't allow any space for it

    textBoxWidth = 0;
    textBoxHeight = 0;
    rightTextMargin = 0;  
    leftTextMargin = 0;
    topTextMargin = 0 ;
    bottomTextMargin = 0;

end

% Set position of bar depending on where we want it
switch(position)
    case 'bottom'
        barX = centreX - barLength / 2;
        barY = height - barWidth - edgeMargin - bottomTextMargin;        
    case 'top'
        barX = centreX - barLength / 2;
        barY = edgeMargin + topTextMargin;
    case 'left'
        barX = edgeMargin + leftTextMargin;
        barY = centreY - barWidth / 2;
    case 'right'
        barX = width - barLength - edgeMargin - rightTextMargin;
        barY = centreY - barWidth / 2;
    case 'top-left'
        barX = edgeMargin + leftTextMargin;
        barY = edgeMargin + topTextMargin;
    case 'top-right'
        barX = width - barLength - edgeMargin - rightTextMargin;
        barY = edgeMargin + topTextMargin;
    case 'bottom-right'
        barX = width - barLength - edgeMargin - rightTextMargin;
        barY = height - barWidth - edgeMargin - bottomTextMargin;
    case 'bottom-left'    
        barX = edgeMargin + leftTextMargin;
        barY = height - barWidth - edgeMargin - bottomTextMargin;
    case 'exact'
       
end


% If we have a text label, find position to write text
if isLabel

    switch(textBoxPosition)
        case 'top'
            textBoxX = barX + barLength / 2 - textBoxWidth / 2;
            textBoxY = barY - textBoxHeight - textBoxPadding;
        case 'bottom'
            textBoxX = barX + barLength / 2 - textBoxWidth / 2;
            textBoxY = barY + barWidth + textBoxPadding;
        case 'left'
            textBoxX = barX - textBoxWidth - textBoxPadding;
            textBoxY = barY + barWidth / 2 - textBoxHeight / 2;
        case 'right'
            textBoxX = barX + barLength + textBoxPadding;
            textBoxY = barY + barWidth / 2 - textBoxHeight / 2;
    end
    
    % This allows for the fact that the text is rendered offset slightly 
    % from the co-ordinates specified when using Matlab's 'text' function 
    textBoxX = textBoxX + textOffsetX;
    textBoxY = textBoxY - textOffsetY;
        
    % Draw a blank image (0) and write the text over it in the correct place
    % so we can work out which pixels we need to copy to the image
    blank(:) = 0;
    imshow(blank);
    t = text(textBoxX, textBoxY, label, 'FontName' , fontName, 'FontSize', fontSize, 'FontWeight', fontWeight, 'FontAngle', fontAngle, 'color', [1 1 1], 'VerticalAlignment', 'top', 'Margin', 1);
    fr = getframe(gca,[1,1,size(image,2),size(image,1)]);
    delete(t);
    c = fr.cdata(:,:,1);
    [textX, textY] = find(c~=0);
    ind = sub2ind(size(blank),textX, textY);
  
end

for iPlane = 1:size(image,3)
   
    % Draw the bar
    image(round(barY):round(barY + barWidth - 1),round(barX):round(barX + barLength -1),iPlane) = colour(iPlane);

end

if isLabel    % If we are adding a text label
    
    % Write the text onto a copy of the image in the correct location
    imshow(image);
    t = text(textBoxX, textBoxY, label,'FontName' , fontName, 'FontSize', fontSize, 'FontWeight', fontWeight, 'FontAngle', fontAngle, 'color', textcolour, 'VerticalAlignment', 'top', 'Margin', 1);
    fr = getframe(gca,[0,0,size(image,2),size(image,1)]);
    delete(t); delete(h)
    
    for iPlane = 1:size(image,3)   % For each colour plane
        
        % Pull out the image data with the text
        c = fr.cdata(:,:,iPlane);
        
        % Pull out the colour plane, copy on the text and then re-insert
        im = image(:,:,iPlane);
        im(ind) = double(c(ind)) /(255 * colourScale);
        image(:,:,iPlane) = im ;
        
    end
    
end
