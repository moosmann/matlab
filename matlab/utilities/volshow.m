function volshow(vol, direction, pause_seconds, fit_to_screen, scroll_by_buttonpress)
% Show slices of a volume using a figure loop
%
% direction: scalar, default: 3. Direction along to scroll.
% pause_seconds: scalar, default: 0.2. Seconds to wait between image
% updates.
% fit_to_screen: transpose image if height of image is larger than its
% width.
% scroll_by_buttonpress: bool, default: 0. Advance loop by pressing
% (almost) any button or clicking in the figure window
%
% Written by Julian Moosmann. First verion: 2016-10-20. Last modification:
% 2016-10-20
%
% volshow(vol, direction, pause_seconds, fit_to_screen, scroll_by_buttonpress)

if nargin < 2
    direction = 3;
end
if nargin < 3
    pause_seconds = 0.1;
end
if nargin < 4
    fit_to_screen = 1;
end
if nargin < 5
    scroll_by_buttonpress = 0;
end

%% Main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

im_shape = size( vol );
im_shape( direction ) = [];
num_sli = size( vol, direction );

% Transpose image if height is larger than width
if fit_to_screen
    if im_shape(1) > im_shape(2)
        fit_to_screen = 1;
        im_shape = im_shape';
    else
        fit_to_screen = 0;
    end
end

% Show loop
show_slice = @(im) imshow( im, [], 'InitialMagnification', 'fit');
fprintf(' Abort with Ctrl-C\n')
for nn = 1:num_sli
    % Update figure name
    fig_str = sprintf(' Slice no %u of %u. Image shape: %u x %u', nn, num_sli, im_shape );
    set( gcf, 'Name', fig_str )
    % Pick slice
    switch direction
        case 1
            im = squeeze( vol(nn, :, :) );            
        case 2
            im = squeeze( vol(:, nn, :) );
        case 3
            im = squeeze( vol(:, :, nn) );
    end
        
    % Transpose?
    show_slice( rot90( FilterHisto(im), fit_to_screen ) )
    
    % Inlay
    x = round( 0.1 * min( size( im )));
    t = text(x, x, sprintf( '%u', nn));
    t.Color = [1 1 1];
    t.FontSize = [16];
    
    drawnow
    shg
    % Interactive scrolling
    if scroll_by_buttonpress
        waitforbuttonpress
    else
        pause( pause_seconds)
    end
end
